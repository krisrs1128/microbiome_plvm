
var width = 800;
var height = 300;

var elem = d3.select("body")
  .append("svg")
  .attrs({
    "width": width,
    "height": height
  });

elem.append("rect")
  .attrs({
    "width": width,
    "height": height,
    "fill": "#F7F7F7"
  });

elem.append("g")
  .attrs({
    "id": "text_box",
    "transform": "translate(" + 0.85 * width + "," + 0.85 * height + ")"
  });

function beta_extract(id) {
  return beta.map(function(x) { return x[id]; });
}

var unique_fill = ["Lachnospiraceae", "Ruminococcaceae", "Bacteroidaceae", "uncultured",
                   "Eubacteria", "Peptostreptococcaceae", "Streptococcaceae", "other"];
var unique_ix = d3.set(beta_extract("ix")).values().map(parseFloat);
var unique_topics = d3.set(beta_extract("topic")).values();

var scales = {
  "fill": d3.scaleOrdinal()
    .domain(unique_fill)
    .range(['#66c2a5','#fc8d62','#8da0cb','#e78ac3','#a6d854','#ffd92f','#e5c494','#b3b3b3']),
  "y": d3.scaleLinear()
    .domain([d3.min(beta_extract("lower")), d3.max(beta_extract("upper"))]),
  "x": d3.scaleLinear()
    .domain(d3.extent(unique_ix))
    .range([0, width]),
  "panels": d3.scaleBand()
    .domain(unique_topics)
    .range([0, 0.85 * height])
};

scales.y.range([scales.panels.step(), 0]);

elem.selectAll("circle")
  .data(beta)
  .enter()
  .append("circle")
  .attrs({
    "class": function(d) {
      var class_text = "circle_" + d.ix + "_" + d.topic;
      return class_text.replace(/\s/g, "_");
    },
    "r": 1,
    "cx": function(d) { return scales.x(d.ix); },
    "cy": function(d) { return scales.panels(d.topic) + scales.y(d.median); },
    "fill": function(d) { return scales.fill(d.fill); }
  });

elem.append("g")
  .attr("id", "error_bars")
  .selectAll("rect")
  .data(beta)
  .enter()
  .append("rect")
  .attrs({
    "class": function(d) {
      var class_text = "error_bar_" + d.ix + "_" + d.topic;
      return class_text.replace(/\s/g, "_");
    },
    "width": 0.2,
    "height": function(d) { return scales.y(d.lower) - scales.y(d.upper); },
    "x": function(d) { return scales.x(d.ix) - 0.1; },
    "y": function(d) { return scales.panels(d.topic) + scales.y(d.upper); },
    "fill": function(d) { return scales.fill(d.fill); }
  });

var voronoi = d3.voronoi()
    .x(function(d) { return scales.x(d.ix); })
    .y(function(d) { return scales.panels(d.topic) + scales.y(d.median); })
    .extent([[0, 0], [width, height]]);

var poly = voronoi(beta)
    .polygons();
poly = poly.filter(function(d) { return !typeof(d) != "undefined"; });

elem.selectAll("path")
	.data(poly)
	.enter()
  .append("path")
  .attr("d", function(d, i) { return "M" + d.join("L") + "Z"; })
  .attrs({
    "class": function(d, i) { return "voronoi" + d.data.ix + "-" + d.data.topic; },
    "fill": "none",
    "pointer-events": "all"
  })
	.on("mouseover", info_over)
	.on("mouseout", info_out);

function info_over(d) {
  d3.select("#text_box")
    .append("text")
    .text(d.data.label)
    .attrs({
      "transform": "translate(10, 20)",
      "font-size": 20
    });

  var select_text = d.data.ix + "_" + d.data.topic;
  select_text = select_text.replace(/\s/g, "_");
  d3.select(".circle_" + select_text)
    .transition()
    .duration(100)
    .attrs({
      "r": 2
    });

  d3.select("#error_bars .error_bar_" + select_text)
    .transition()
    .duration(100)
    .attrs({
      "width": 1.8,
      "x": function(d) { return scales.x(d.ix) - 0.9; }
    });
}

function info_out(d) {
  d3.select("#text_box")
    .selectAll("text")
    .remove();

  var select_text = d.data.ix + "_" + d.data.topic;
  select_text = select_text.replace(/\s/g, "_");
  d3.select(".circle_" + select_text)
    .transition()
    .duration(100)
    .attrs({
      "r": 1
    });

  d3.select("#error_bars .error_bar_" + select_text)
    .transition()
    .duration(100)
    .attrs({
      "width": 0.2,
      "x": function(d) { return scales.x(d.ix) - 0.1; }
    });
}
