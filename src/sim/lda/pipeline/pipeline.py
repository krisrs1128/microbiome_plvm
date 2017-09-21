"""
Pipeline for Evaluating LDA inferences on Simulated Data

This is a luigi pipeline [http://luigi.readthedocs.io/] pipeline for simulating
data from the LDA model and evaluating inference based on three procedures --
Gibbs Sampling, Variational Bayes, and the Bootstrap. Running

python3 pipeline.py LDAExperiment --local-scheduler --workers=5

will start three processes in parallel, writing simulation output to the
directory specified by output_dir in luigi.cfg. The simulation parameters must
be specified in a configuration file, also linked through the luigi.cfg.

Note however that this script does none of the post-simulation visualization.
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
logger = logging.getLogger("lda.pipeline")

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
        "".join([self.K_fit, self.alpha_fit, self.gamma_fit, self.D,
                 self.N, self.V, self.K, self.alpha0, self.gamma0,
                 self.conf.get("expers", "n_samples")])
    )


def run_and_check(cmds):
    run_cmd = [str(s) for s in cmds]
    status = subprocess.call(run_cmd)
    if status is not 0:
        raise ValueError("Bash commands failed")

###############################################################################
# Core pipeline classes
###############################################################################

class LDAExperiment(luigi.WrapperTask):
    """Wrapper Experiment Task

    Run a complete LDA simulation experiment. This wraps all the gibbs,
    bootstrap, and variational bayes inference procedures. It loops over the
    experiments configuration file and requires the variational bayes and
    bootstrapping output for combination of simulation parameter settings.

    Args: None

    Attributes:
        conf (configuration): A luigi configuration object, created by parsing
        ./luigi.cfg. This provides the link to high level experiment
        parameters.
    """
    conf = configuration.get_config()

    def requires(self):
        """
        get the experiments configuration file
        """
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
                str(v["K"]), str(v["alpha0"]), str(v["gamma0"]),
                str(v["D"]), str(v["N"]), str(v["V"]), str(v["K"]),
                str(v["alpha0"]), str(v["gamma0"])
            ]

            for i, _ in enumerate(batch_endpoints[:-1]):
                boot_params = [batch_endpoints[i], batch_endpoints[i + 1] - 1] + \
                              data_params
                tasks.append(LDABoot(*boot_params))

            gibbs_params = ["gibbs"] + data_params
            tasks.append(LDAFit(*gibbs_params))

        return tasks


class LDABoot(luigi.Task):
    """
    Parametric Bootstrap inference for LDA

    This task generates parametric bootstrap samples from a fitted LDA model.
    If the original VB fit has not been generated, this will launch that task.

    Arguments:
      start_ix (int): The start index for the current bootstrapping iterations
      end_ix (int): The end index for the current bootstrapping iterations
      K_fit (int): How many topics will we tell LDA to use when fitting?
      alpha_fit (float): What is the theta parameter prior we should use
        across all coordinates, for fitting?
      gamma_fit (float): What is the beta parameter prior we should use
        across all coordinates, for fitting?
      D (int): How many samples are there in this experiment?
      N (int): How many words are there in each sample?
      V (int): How many terms are there across samples?
      K (int): How many topics are there?
      alpha0 (float): What is the true theta parameter prior used in
        generating data?
      gamma0 (float): What is the true beta parameter prior used in generating
        data?

    Attributes:
        See arguments
    """
    start_ix = luigi.Parameter()
    end_ix = luigi.Parameter()
    K_fit = luigi.Parameter()
    alpha_fit = luigi.Parameter()
    gamma_fit = luigi.Parameter()
    D = luigi.Parameter()
    N = luigi.Parameter()
    V = luigi.Parameter()
    K = luigi.Parameter()
    alpha0 = luigi.Parameter()
    gamma0 = luigi.Parameter()

    conf = configuration.get_config()

    def requires(self):
        return LDAFit(
            "vb", self.K_fit, self.alpha_fit, self.gamma_fit, self.D,
            self.N, self.V, self.K, self.alpha0, self.gamma0
        )

    def run(self):
        input_path = self.input().open("r").name

        run_cmd = [
            "Rscript", self.conf.get("expers", "boot_script"),
            os.path.join(self.conf.get("expers", "output_dir"), "bootstraps"),
            self.start_ix, self.end_ix, fit_id(self), self.conf.get("expers",
            "n_samples"), input_path, self.N, self.alpha0, self.gamma0,
        ]
        print(run_cmd)
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

        theta_paths = [luigi.LocalTarget(os.path.join(output_dir, "theta-" + s)) for s in output_base]
        beta_paths = [luigi.LocalTarget(os.path.join(output_dir, "beta-" + s)) for s in output_base]

        return theta_paths + beta_paths


class LDAFit(luigi.Task):
    """
    Fit an LDA model on the simulated data .

    This task fits an LDA (gibbs or VB) model to simulated LDA data, or
    triggers the data simulation if it is not available.

    Arguments:
      fit_method (str): Should we use variational bayes or gibbs sampling?
      K_fit (int): How many topics will we tell LDA to use when fitting?
      alpha_fit (float): What is the theta parameter prior we should use
        across all coordinates, for fitting?
      gamma_fit (float): What is the beta parameter prior we should use
        across all coordinates, for fitting?
      D (int): How many samples are there in this experiment?
      N (int): How many words are there in each sample?
      V (int): How many terms are there across samples?
      K (int): How many topics are there?
      alpha0 (float): What is the true theta parameter prior used in
        generating data?
      gamma0 (float): What is the true beta parameter prior used in generating
        data?

    Attributes:
        See arguments
    """
    fit_method = luigi.Parameter()
    K_fit = luigi.Parameter()
    alpha_fit = luigi.Parameter()
    gamma_fit = luigi.Parameter()
    D = luigi.Parameter()
    N = luigi.Parameter()
    V = luigi.Parameter()
    K = luigi.Parameter()
    alpha0 = luigi.Parameter()
    gamma0 = luigi.Parameter()

    conf = configuration.get_config()

    def requires(self):
        return LDAData(self.D, self.N, self.V, self.K, self.alpha0, self.gamma0)

    def run(self):
        data_path = self.input().open("r").name

        run_cmd = [
            "Rscript", self.conf.get("expers", "fit_script"),
            self.conf.get("expers", "output_dir"), fit_id(self), data_path,
            self.fit_method, self.conf.get("expers", "n_samples"), self.K_fit,
            self.alpha_fit, self.gamma_fit
        ]
        run_and_check(run_cmd)

    def output(self):
        output_dir = self.conf.get("expers", "output_dir")
        output_file = self.fit_method + "-" + fit_id(self) + ".RData"
        return luigi.LocalTarget(os.path.join(output_dir, output_file))


class LDAData(luigi.Task):
    """
    Simulate Data from Theta and Beta parameters

    This task generates data according to the LDA model, assuming the
    parameters exist, or otherwise triggers the task for simulating data
    parameters.

    Arguments:
      D (int): How many samples are there in this experiment?
      N (int): How many words are there in each sample?
      V (int): How many terms are there across samples?
      K (int): How many topics are there?
      alpha0 (float): What is the true theta parameter prior used in
        generating data?
      gamma0 (float): What is the true beta parameter prior used in generating
        data?

    Attributes:
      see arguments
    """
    D = luigi.Parameter()
    N = luigi.Parameter()
    V = luigi.Parameter()
    K = luigi.Parameter()
    alpha0 = luigi.Parameter()
    gamma0 = luigi.Parameter()

    conf = configuration.get_config()

    def requires(self):
        return LDAParams(self.D, self.V, self.K, self.alpha0, self.gamma0)

    def run(self):
        beta_path = self.input()[0].open("r").name
        theta_path = self.input()[1].open("r").name
        gen_id = hash_string("".join([self.D, self.N, self.V, self.K, self.alpha0, self.gamma0]))

        run_cmd = [
            "Rscript",
            self.conf.get("expers", "sim_script"),
            self.conf.get("expers", "output_dir"),
            gen_id, self.N, beta_path, theta_path
        ]
        run_and_check(run_cmd)

    def output(self):
        gen_id = hash_string("".join([self.D, self.N, self.V, self.K, self.alpha0, self.gamma0]))
        output_dir = self.conf.get("expers", "output_dir")
        return luigi.LocalTarget(os.path.join(output_dir, "n-" + gen_id + ".feather"))


class LDAParams(luigi.ExternalTask):
    """
    Simulate parameters for an LDA model

    This generates the beta and theta parameters that are the basis for the rest of the simulation.

    Arguments:
      D (int): How many samples are there in this experiment?
      V (int): How many terms are there across samples?
      K (int): How many topics are there?
      alpha0 (float): What is the true theta parameter prior used in
        generating data?
      gamma0 (float): What is the true beta parameter prior used in generating
        data?
    """
    D = luigi.Parameter()
    V = luigi.Parameter()
    K = luigi.Parameter()
    alpha0 = luigi.Parameter()
    gamma0 = luigi.Parameter()

    conf = configuration.get_config()

    def run(self):
        gen_id = hash_string("".join([self.D, self.V, self.K, self.alpha0, self.gamma0]))
        run_cmd = [
            "Rscript", self.conf.get("expers", "param_script"),
            self.conf.get("expers", "output_dir"), gen_id, self.D, self.V,
            self.K, self.alpha0, self.gamma0
        ]
        run_and_check(run_cmd)

    def output(self):
        gen_id = hash_string("".join([self.D, self.V, self.K, self.alpha0, self.gamma0]))
        output_dir = self.conf.get("expers", "output_dir")
        return [
            luigi.LocalTarget(os.path.join(output_dir, "beta-" + gen_id + ".feather")),
            luigi.LocalTarget(os.path.join(output_dir, "theta-" + gen_id + ".feather"))
        ]


if __name__ == "__main__":
    luigi.run()
