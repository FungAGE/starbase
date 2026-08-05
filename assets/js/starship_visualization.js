/**
 * Starship gene-feature map - D3 v6, ported from MAS4starships'
 * viz-upgrade branch (starship_visualization.js, D3 v7 rewrite).
 *
 * Adapted for Dash: no DOMContentLoaded/json_script -- exposes
 * window.renderStarshipViz(containerId, data) instead, called from a
 * clientside_callback whenever the reviewed annotation's ship changes.
 * data: {starship_length, features: [{id,start,stop,strand,type,
 *   annotation_id,flag_label,annotation,public_note}], current_annotation_id}
 */

(function () {
  'use strict';

  const FEATURE_DENSITY_THRESHOLD = 800;
  const NUM_DENSITY_BINS = 80;

  function range(start, stop, step = 1) {
    return Array(Math.ceil((stop - start) / step))
      .fill(start)
      .map((x, y) => x + y * step);
  }

  // Keyed by FLAG_LABELS values from src/database/curation_constants.py
  const FLAG_COLORS = {
    GREEN: 'green',
    YELLOW: '#f1c40f',
    RED: 'red',
    REVIEW_NAME: '#C21460',
    N_A: 'grey',
    ORANGE: 'orange',
    UNANNOTATED: 'grey',
    CURRENT: 'fuchsia',
  };

  function getFeatureColor(d, currentAnnotationId) {
    if (currentAnnotationId && d.annotation_id === currentAnnotationId) {
      return FLAG_COLORS.CURRENT;
    }
    if (d.annotation === 'hypothetical protein') return 'black';
    return FLAG_COLORS[d.flag_label] || '#888';
  }

  function bottom_strand(previous, start, end) {
    if (previous['end'] !== null) {
      if (start < previous.end) {
        previous['end'] = end;
        return previous['x'] === '-25' ? (previous['x'] = '-30') : (previous['x'] = '-25');
      }
      previous['end'] = end;
      previous['x'] = '-25';
      return '-25';
    }
    previous['x'] = '-25';
    return '-25';
  }

  function top_strand(previous, start, end) {
    if (previous['end'] !== null) {
      if (end < previous.end) {
        previous['end'] = start;
        return previous['x'] === '-50' ? (previous['x'] = '-55') : (previous['x'] = '-50');
      }
      previous['end'] = start;
      previous['x'] = '-50';
      return '-50';
    }
    previous['x'] = '-50';
    return '-50';
  }

  function render(containerId, data) {
    const container = document.getElementById(containerId);
    if (!container || !data) return;

    container.innerHTML = '';

    const features = data.features || [];
    const starshipLength = data.starship_length || 1;
    const currentAnnotationId = data.current_annotation_id || null;

    if (!features.length) {
      container.innerHTML =
        '<div style="padding: 1rem; color: var(--mantine-color-gray-6);">No gene features on this ship yet.</div>';
      return;
    }

    const margin = { top: 50, right: 20, bottom: 50, left: 60 };
    const containerWidth = container.getBoundingClientRect().width || 800;
    const containerHeight = 260;
    const width = Math.max(100, containerWidth * 0.95 - margin.left - margin.right);
    const height = Math.max(150, containerHeight - margin.top - margin.bottom);

    const svg = d3
      .select('#' + containerId)
      .append('svg')
      .attr('width', '100%')
      .attr('height', containerHeight)
      .attr('class', 'starship-viz-svg');

    const xScale = d3.scaleLinear().domain([0, starshipLength]).range([0, width]);

    const scaleVertical = d3.scaleLinear().domain([2, 0]).range([0, 50]);
    const yAxis = d3.axisLeft(scaleVertical).tickFormat((d, i) => {
      if (i === 0) return '0';
      return i % 2 === 0 ? '-' : '+';
    });
    const xAxis = d3.axisBottom().scale(xScale);

    const zoom = d3
      .zoom()
      .scaleExtent([1, 100])
      .on('zoom', (event) => zoomed(event));

    const chartGroup = svg.append('g').attr('transform', `translate(${margin.left},${margin.top})`);

    const starshipVis = chartGroup
      .append('g')
      .attr('id', containerId + '-vis')
      .attr('transform', `translate(0, ${height - 80})`)
      .call(xAxis);

    chartGroup
      .append('text')
      .attr('transform', `translate(${width / 2}, ${height - 50})`)
      .style('text-anchor', 'middle')
      .text('Starship');

    chartGroup
      .append('g')
      .attr('class', 'yAxis')
      .attr('transform', 'translate(0, 80)')
      .call(yAxis.tickValues(range(0, 3, 1)));

    chartGroup
      .append('text')
      .attr('transform', 'rotate(-90)')
      .attr('y', 0)
      .attr('x', -(height / 2))
      .attr('dy', '1em')
      .style('text-anchor', 'middle')
      .text('Strand');

    const zoomControls = chartGroup
      .append('g')
      .attr('class', 'starship-zoom-controls')
      .attr('transform', `translate(${width - 80}, 5)`)
      .style('cursor', 'pointer');

    zoomControls
      .append('rect')
      .attr('width', 75)
      .attr('height', 28)
      .attr('rx', 4)
      .attr('fill', '#f5f5f5')
      .attr('stroke', '#999');

    zoomControls
      .append('text')
      .attr('x', 37)
      .attr('y', 18)
      .style('text-anchor', 'middle')
      .style('font-size', '12px')
      .style('pointer-events', 'none')
      .text('Reset zoom');

    zoomControls.on('click', () => {
      svg.transition().duration(300).call(zoom.transform, d3.zoomIdentity);
    });

    let tooltip = d3.select('body').select('.starship-feature-tooltip');
    if (tooltip.empty()) {
      tooltip = d3
        .select('body')
        .append('div')
        .attr('class', 'starship-feature-tooltip')
        .style('position', 'absolute')
        .style('visibility', 'hidden')
        .style('background', 'rgba(0,0,0,0.85)')
        .style('color', '#fff')
        .style('padding', '8px 12px')
        .style('border-radius', '4px')
        .style('font-size', '12px')
        .style('max-width', '320px')
        .style('z-index', '1000')
        .style('pointer-events', 'none')
        .style('white-space', 'pre-line');
    }

    function showTooltip(event, d) {
      const text = [
        `Feature: ${d.start}–${d.stop} (${d.strand || 'N/A'})`,
        `Type: ${d.type || 'N/A'} | Flag: ${d.flag_label || 'N/A'}`,
        `Annotation: ${d.annotation || 'N/A'}`,
        d.public_note ? `Note: ${d.public_note}` : null,
      ]
        .filter(Boolean)
        .join('\n');
      tooltip.style('visibility', 'visible').text(text);
    }

    function moveTooltip(event) {
      tooltip.style('left', event.pageX + 12 + 'px').style('top', event.pageY + 12 + 'px');
    }

    function hideTooltip() {
      tooltip.style('visibility', 'hidden');
    }

    function drawFeatures(scale, visGroup) {
      const domain = scale.domain();
      const viewStart = Math.max(0, Math.floor(domain[0]));
      const viewEnd = Math.min(starshipLength, Math.ceil(domain[1]));
      const inView = features.filter(
        (f) =>
          (f.start <= viewEnd && f.stop >= viewStart) ||
          (f.start > f.stop && (f.start <= viewEnd || f.stop >= viewStart))
      );

      visGroup.selectAll('.feature-rect').remove();
      visGroup.selectAll('.density-bar').remove();

      if (inView.length > FEATURE_DENSITY_THRESHOLD) {
        const binWidth = (viewEnd - viewStart) / NUM_DENSITY_BINS || 1;
        const bins = Array(NUM_DENSITY_BINS)
          .fill(0)
          .map((_, i) => ({ start: viewStart + i * binWidth, end: viewStart + (i + 1) * binWidth, count: 0 }));
        inView.forEach((f) => {
          const start = Math.min(f.start, f.stop);
          const end = Math.max(f.start, f.stop);
          const binStart = Math.floor((start - viewStart) / binWidth);
          const binEnd = Math.min(NUM_DENSITY_BINS - 1, Math.ceil((end - viewStart) / binWidth));
          for (let b = Math.max(0, binStart); b <= binEnd; b++) bins[b].count++;
        });

        const maxCount = d3.max(bins, (b) => b.count) || 1;
        const barHeight = 20;
        visGroup
          .selectAll('.density-bar')
          .data(bins.filter((b) => b.count > 0))
          .join('rect')
          .attr('class', 'density-bar')
          .attr('x', (d) => scale(d.start))
          .attr('y', (d) => 50 - (d.count / maxCount) * barHeight)
          .attr('width', (d) => Math.max(1, scale(d.end) - scale(d.start)))
          .attr('height', (d) => Math.max(1, (d.count / maxCount) * barHeight))
          .attr('fill', '#69b3a2')
          .attr('opacity', 0.7)
          .append('title')
          .text((d) => `${d.count} features in ${Math.round(d.end - d.start)} bp`);
      } else {
        const previousb = {};
        const previoust = {};
        visGroup
          .selectAll('.feature-rect')
          .data(inView)
          .join('rect')
          .attr('class', 'feature-rect')
          .attr('x', (d) => (d.strand === '+' ? scale(d.start) : scale(d.stop)))
          .attr('y', (d) =>
            d.strand === '+'
              ? parseFloat(bottom_strand(previousb, d.start, d.stop))
              : parseFloat(top_strand(previoust, d.start, d.stop))
          )
          .attr('width', (d) => Math.abs(scale(Math.max(d.start, d.stop)) - scale(Math.min(d.start, d.stop))))
          .attr('height', 5)
          .style('fill', (d) => getFeatureColor(d, currentAnnotationId))
          .style('opacity', 0.9)
          .style('cursor', 'pointer')
          .on('mouseover', (event, d) => {
            showTooltip(event, d);
            d3.select(event.currentTarget).style('opacity', 1);
          })
          .on('mousemove', moveTooltip)
          .on('mouseout', (event) => {
            hideTooltip();
            d3.select(event.currentTarget).style('opacity', 0.9);
          });
      }
    }

    const visGroup = starshipVis.append('g').attr('transform', 'translate(0, -100)');
    drawFeatures(xScale, visGroup);

    function zoomed(event) {
      const newScale = event.transform.rescaleX(xScale);
      starshipVis.transition().duration(0).call(xAxis.scale(newScale));
      drawFeatures(newScale, visGroup);
    }

    svg.call(zoom);

    const legendData = [
      { color: 'black', label: 'Hypothetical protein' },
      { color: 'green', label: 'Green' },
      { color: '#f1c40f', label: 'Yellow' },
      { color: 'red', label: 'Red' },
      { color: '#C21460', label: 'Review name' },
      { color: 'orange', label: 'Orange' },
      { color: 'grey', label: 'N/A, Unannotated' },
      { color: 'fuchsia', label: 'Current annotation' },
    ];

    const legend = chartGroup
      .append('g')
      .attr('class', 'starship-legend')
      .attr('transform', `translate(0, ${height - 70})`);

    legend
      .selectAll('.legend-item')
      .data(legendData)
      .join('g')
      .attr('class', 'legend-item')
      .attr('transform', (_, i) => `translate(${Math.floor(i / 4) * 220}, ${(i % 4) * 18})`)
      .each(function (d) {
        const g = d3.select(this);
        g.append('circle').attr('r', 5).attr('cx', 0).attr('cy', -5).style('fill', d.color);
        g.append('text').attr('x', 10).attr('y', 0).style('font-size', '11px').text(d.label);
      });
  }

  window.renderStarshipViz = render;
})();
