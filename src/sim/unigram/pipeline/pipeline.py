"""
Pipeline for Evaluating Unigram inferences on Simulated Data
"""

import logging
import logging.config

import luigi
from luigi import configuration
import subprocess
import os
import hashlib
import json

logging_conf = configuration.get_config().get("core", "logging_conf_file")
logging.config.fileConfig(logging_conf)
logger = logging.getLogger("unigram.pipeline")

###############################################################################
# Minor utilities used throughout the pipeline
###############################################################################


def hash_string(string, max_chars=32):
    """
    Return the sha1 hash of a string
    """
    hash_obj = hashlib.sha1(string.encode())
    return hash_obj.hexdigest()[:max_chars]


def fit_id(self):
    return hash_string(
        "".join([self.conf.get("expers", "n_samples"), self.a0, self.b0,
                 self.D, self.N, self.V, self.sigma0])
    )


def run_and_check(cmds):
    run_cmd = [str(s) for s in cmds]
    status = subprocess.call(run_cmd)
    if status is not 0:
        raise ValueError("Bash commands failed")


###############################################################################
# Core pipeline classes
###############################################################################
class UnigramExperiment(luigi.WrapperTask):
    """Wrapper Task for Full Unigram Experiment

    Run the complete Unigram simulation experiment. It will send off the gibbs,
    bootstrap, and variational bayes runs that are later visualized in the
    simulations section.

    Attributes:
        conf (configuration): A luigi configuration object, created by parsing
        ./luigi.cfg. This provides the link to high level experiment
        parameters.
    """
    conf = configuration.get_config()

    def requires(self):
        with open(self.conf.get("expers", "master")) as f:
            experiment = json.load(f)

        n_batches = int(self.conf.get("expers", "n_batches"))
        n_samples = int(self.conf.get("expers", "n_samples"))
        batches = (list(range(n_batches)) * (int(n_samples / n_batches) + 1))[:n_samples]
        batches.sort()
        batch_endpoints = [0] + [
            i for i in list(range(n_samples - 1))
            if batches[i] != batches[i + 1]
        ]

        tasks = []
        for (k, v) in enumerate(experiment):
            data_params = [
                str(v["a0"]), str(v["b0"]), str(v["D"]), str(v["N"]),
                str(v["V"]), str(v["sigma0"])
            ]

            for i, _ in enumerate(batch_endpoints[:-1]):
                boot_params = [batch_endpoints[i], batch_endpoints[i + 1] - 1] + \
                              data_params
                tasks.append(UnigramBoot(*boot_params))

            gibbs_params = ["gibbs"] + data_params
            tasks.append(UnigramFit(*gibbs_params))

        return tasks


class UnigramBoot(luigi.Task):
    """
    Parametric Boostrap inference for Dynamic Unigrams

    Arguments:
      start_ix (int): The start index for the current bootstrapping iterations
      end_ix (int): The end index for hte current bootstrapping interations
      a0 (float): What is the hyperparameter a0 we should use when fitting?
      b0 (float): What is the hyperparameter b0 we should use when fitting?
      D (int): How many samples are there in this experiment?
      N (int): How many words are there in each sample?
      V (int): How many terms are there across samples?
    """
    start_ix = luigi.Parameter()
    end_ix = luigi.Parameter()
    a0 = luigi.Parameter()
    b0 = luigi.Parameter()
    D = luigi.Parameter()
    N = luigi.Parameter()
    V = luigi.Parameter()
    sigma0 = luigi.Parameter()

    conf = configuration.get_config()

    def requires(self):
        return UnigramFit(
            "vb", self.a0, self.b0, self.D, self.N, self.V, self.sigma0
        )

    def run(self):
        input_path = self.input().open("r").name
        run_cmd = [
            "Rscript", self.conf.get("expers", "boot_script"),
            os.path.join(self.conf.get("expers", "output_dir"), "bootstraps"),
            self.start_ix, self.end_ix, fit_id(self), input_path,
            self.conf.get("expers", "stan_path"), self.conf.get("expers",
            "n_samples"), self.N, self.a0, self.b0,
        ]
        run_and_check(run_cmd)

    def output(self):
        output_dir = os.path.join(
            self.conf.get("expers", "output_dir"),
            "bootstraps"
        )
        output_base = [
            "boot-" + fit_id(self) + str(i) + ".feather"
            for i in range(int(self.start_ix), int(self.end_ix))
        ]

        return [luigi.LocalTarget(os.path.join(output_dir, "mu-" + s))
                for s in output_base]


class UnigramFit(luigi.Task):
    """
    Fit a Unigram model on simulated data

    Arguments:
      fit_method (str): Should we use variational bayes or gibbs sampling?
      a0 (float): What is the hyperparameter a0 we should use when fitting?
      b0 (float): What is the hyperparameter b0 we should use when fitting?
      D (int): How many samples are there in this experiment?
      N (int): How many words are there in each sample?
      V (int): How many terms are there across samples?
      sigma0 (float): What is the true sigma random walk size parameter used in
      generating the data?
    """

    fit_method = luigi.Parameter()
    a0 = luigi.Parameter()
    b0 = luigi.Parameter()
    D = luigi.Parameter()
    N = luigi.Parameter()
    V = luigi.Parameter()
    sigma0 = luigi.Parameter()

    conf = configuration.get_config()

    def requires(self):
        return UnigramData(self.D, self.N, self.V, self.sigma0)

    def run(self):
        run_cmd = [
            "Rscript",
            self.conf.get("expers", "fit_script"),
            self.conf.get("expers", "output_dir"),
            self.fit_method,
            self.conf.get("expers", "stan_path"),
            fit_id(self),
            self.input().open("r").name,
            self.conf.get("expers", "n_samples"),
            self.a0,
            self.b0
        ]
        run_and_check(run_cmd)

    def output(self):
        output_dir = self.conf.get("expers", "output_dir")
        output_file = self.fit_method + "-" + fit_id(self) + ".RData"
        return luigi.LocalTarget(os.path.join(output_dir, output_file))


class UnigramData(luigi.Task):
    """
    Simulate data according to a Unigram model

    Arguments:
      D (int): How many samples are there in this experiment?
      N (int): How many words are there in each sample?
      V (int): How many terms are there across samples?
      sigma0 (float): What is the true sigma random walk size parameter used in
      generating the data?
    """
    D = luigi.Parameter()
    N = luigi.Parameter()
    V = luigi.Parameter()
    sigma0 = luigi.Parameter()

    conf = configuration.get_config()

    def requires(self):
        return UnigramParams(self.D, self.V, self.sigma0)

    def run(self):
        mu_path = self.input()[0].open("r").name
        gen_id = hash_string("".join([self.D, self.N, self.V, self.sigma0]))
        run_cmd = [
            "Rscript",
            self.conf.get("expers", "sim_script"),
            self.conf.get("expers", "output_dir"),
            gen_id,
            self.N,
            mu_path
        ]
        run_and_check(run_cmd)

    def output(self):
        gen_id = hash_string("".join([self.D, self.N, self.V, self.sigma0]))
        output_dir = self.conf.get("expers", "output_dir")
        return luigi.LocalTarget(os.path.join(output_dir, "x-" + gen_id + ".feather"))


class UnigramParams(luigi.ExternalTask):
    """
    Simulate parameters for a Unigram model

    This generates the parameters mu[t] for a particular instance of the
    dynamic unigram model.

    Arguments:
      D (int): How many samples are there in this experiment?
      V (int): How many terms are there across samples?
      sigma0 (float): What is the true sigma random walk size parameter used in
      generating the data?
    """
    D = luigi.Parameter()
    V = luigi.Parameter()
    sigma0 = luigi.Parameter()

    conf = configuration.get_config()

    def run(self):
        print("test")
        gen_id = hash_string("".join([self.D, self.V, self.sigma0]))
        print(gen_id)
        run_cmd = [
            "Rscript", self.conf.get("expers", "param_script"),
            self.conf.get("expers", "output_dir"), gen_id, self.D, self.V,
            self.sigma0
        ]
        run_and_check(run_cmd)

    def output(self):
        gen_id = hash_string("".join([self.D, self.V, self.sigma0]))
        output_dir = self.conf.get("expers", "output_dir")
        return [
            luigi.LocalTarget(os.path.join(output_dir, "mu-" + gen_id + ".feather"))
        ]


if __name__ == "__main__":
    luigi.run()
