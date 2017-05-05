
var width = 600;
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
    "transform": "translate(" + 0.65 * width + "," + 0.85 * height + ")",
  });

var unique_fill = d3.set(beta.map(function(x) { return x.fill; })).values();
var unique_ix = d3.set(beta.map(function(x) { return x.ix; })).values();
var unique_topics = d3.set(beta.map(function(x) { return x.topic; })).values();

var scales = {
  "fill": d3.scaleOrdinal()
    .domain(unique_fill)
    .range(['#66c2a5','#fc8d62','#8da0cb','#e78ac3']),
  "y": d3.scaleLinear()
    .domain(d3.extent(beta.map(function(x) { return x.median; }))),
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
      var class_text = "circle-" + d.ix + "-" + d.topic;
      return class_text.replace(/\s/g, "-");
    },
    "r": 1,
    "cx": function(d) { return scales.x(d.ix); },
    "cy": function(d) { return scales.panels(d.topic) + scales.y(d.median); },
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
    .text(d.data.ix)
    .attrs({
      "transform": "translate(10, 20)",
      "font-size": 20
    });

  var select_text = ".circle-" + d.data.ix + "-" + d.data.topic;
  select_text = select_text.replace(/\s/g, "-");
  d3.select(select_text)
    .transition()
    .duration(100)
    .attrs({
      "r": 4
    });
}

function info_out(d) {
  d3.select("#text_box")
    .selectAll("text")
    .remove();

  var select_text = ".circle-" + d.data.ix + "-" + d.data.topic;
  select_text = select_text.replace(/\s/g, "-");
  d3.select(select_text)
    .transition()
    .duration(100)
    .attrs({
      "r": 1
    });
}
