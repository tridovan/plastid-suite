/* The page of the local interface. It draws what the server reads from the run directories and
   computes nothing: every number on it is a cell of a table a run wrote. */
'use strict';

const SVG = 'http://www.w3.org/2000/svg';
const STATUS = {
  ok: ['Finished', 'ok'], failed: ['Failed', 'bad'], refused: ['Refused', 'warn'],
  tool_error: ['Tool problem', 'tool'], conflict: ['Conflict', 'conflict'],
  cancelled: ['Cancelled', 'quiet'], interrupted: ['Interrupted', 'quiet'],
  running: ['Running', 'run'], partial: ['Partly finished', 'warn'],
  incomplete: ['Incomplete', 'warn'], unknown: ['Unknown', 'quiet'],
};
const ROUTE = [
  ['assemble', 'Assemble', 'reads to a sequence'],
  ['annotate', 'Annotate', 'a sequence to genes'],
  ['records', 'Records', 'the IR gate'],
  ['analyse', 'Analyse', 'choose analyses'],
  ['runs', 'Results', 'tables, figures, reports'],
  ['export', 'Export', 'package, submission'],
];
const MAPS_SHOWN = 6;  // a run of 34 genomes is not read as 34 strips
// The order a person reads analyses in: the genome first, its parts next, comparisons by question.
const ORDER = ['structure', 'genome_stats', 'gene_table', 'introns', 'orf', 'codon_usage', 'amino_acid',
  'ssr', 'repeats', 'rna_editing', 'submission', 'diversity', 'substitutions', 'barcodes', 'kaks', 'popgen',
  'synteny', 'stats', 'phylogeny', 'gene_names', 'gene_normalize'];
// The cells a run's page leads with: [table file, column, label, unit or note].
const HIGHLIGHTS = [
  ['genome_stats', 'Length (bp)', 'Length', 'bp, whole plastome'],
  ['genome_stats', 'GC (%)', 'GC', '%, whole plastome'],
  ['genome_stats', 'Genes', 'Genes', 'CDS + tRNA + rRNA, every copy'],
  ['ssr_summary', 'Total SSRs', 'SSR', 'both IR copies included'],
  ['repeat_summary', 'Matches', 'Exact repeats', 'exact matches only'],
  ['rna_editing_summary', 'Predicted sites', 'Editing sites', 'predicted, not observed'],
  ['diversity_summary', 'π (whole alignment)', 'π', 'whole alignment, gaps dropped'],
  ['substitutions_summary', 'Substitutions', 'Substitutions', 'see the scope below'],
  ['substitutions_summary', 'Ts/Tv', 'Ts/Tv', 'see the scope below'],
  ['barcodes_summary', 'Hypervariable regions', 'Hypervariable', 'regions, mean + 1 sd'],
  ['kaks_pooled', 'Pooled Ka/Ks', 'Ka/Ks, pooled', 'over the genes scored'],
  ['synteny_genomes', 'Coverage (%)', 'Coverage', '% of the reference'],
];

const app = document.getElementById('app');
const session = { about: null, modules: null, tools: null, chosen: new Set(), viaSequence: new Set(), form: null };
let timers = [];

/* Small tools -------------------------------------------------------------------------------- */

function el(tag, attrs, ...children) {
  const node = document.createElement(tag);
  for (const [key, value] of Object.entries(attrs || {})) {
    if (value === null || value === undefined || value === false) continue;
    if (key.startsWith('on')) node.addEventListener(key.slice(2), value);
    else if (key === 'class') node.className = value;
    else if (key === 'value') node.value = value;
    else if (key === 'checked') node.checked = true;
    else node.setAttribute(key, value === true ? '' : value);
  }
  for (const child of children.flat(Infinity)) {
    if (child === null || child === undefined || child === false) continue;
    node.append(child.nodeType ? child : document.createTextNode(String(child)));
  }
  return node;
}

function svg(tag, attrs, ...children) {
  const node = document.createElementNS(SVG, tag);
  for (const [key, value] of Object.entries(attrs || {})) node.setAttribute(key, value);
  node.append(...children);
  return node;
}

async function call(url, options) {
  const response = await fetch(url, options);
  const payload = await response.json().catch(() => ({ error: `the server answered ${response.status}` }));
  if (!response.ok) throw new Error(payload.error || `the server answered ${response.status}`);
  return payload;
}
const get = (path, query) => call(path + (query ? `?${new URLSearchParams(query)}` : ''));
const post = (path, body) => call(path, {
  method: 'POST', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify(body || {}),
});

const every = (ms, work) => { timers.push(setInterval(work, ms)); };
const number = (value) => (value === null || value === undefined ? '' : Number(value).toLocaleString('en-US'));
const statusOf = (status) => STATUS[status] || [status, 'quiet'];
const pill = (status) => el('span', { class: `pill ${statusOf(status)[1]}` }, statusOf(status)[0]);
const word = (status) => el('span', { class: `word-${statusOf(status)[1]}` }, statusOf(status)[0]);
const species = (text) => el('em', { class: 'species' }, text);
const when = (stamp) => (stamp ? stamp.replace('T', ' ').slice(0, 16) : '');

// An answer that arrives late must not paint over the page the person has since moved to.
function page(active, ...content) {
  if (active && active !== session.active) return;
  for (const timer of timers) clearInterval(timer);
  timers = [];
  app.replaceChildren(rail(active), el('main', null, ...content));
  window.scrollTo(0, 0);
}

function fail(error) {
  return el('div', { class: 'note-warn' }, String(error.message || error));
}

/* The rail -------------------------------------------------------------------------------------- */

function ring() {
  const dash = [['#3FBF93', '40.78 34.62', 0], ['#7BA0F5', '12.86 62.54', -40.78],
    ['#F4B04A', '8.91 66.49', -53.64], ['#7BA0F5', '12.86 62.54', -62.55]];
  const group = svg('g', { transform: 'rotate(-90 16 16)', fill: 'none', 'stroke-width': 5 });
  for (const [colour, pattern, offset] of dash) {
    group.append(svg('circle', { cx: 16, cy: 16, r: 12, stroke: colour, 'stroke-dasharray': pattern, 'stroke-dashoffset': offset }));
  }
  return svg('svg', { width: 34, height: 34, viewBox: '0 0 32 32', 'aria-hidden': 'true' }, group);
}

function rail(active) {
  const about = session.about || {};
  const tools = session.tools ? session.tools.tools : null;
  const found = tools ? tools.filter((tool) => tool.status === 'ok').length : null;
  return el('nav', { class: 'rail', 'aria-label': 'Main' },
    el('a', { class: 'brand', href: '#/start' }, ring(),
      el('span', null, el('b', null, 'PlastidToolSuite'), el('small', null, `version ${about.version || ''}`))),
    el('a', { class: `rail-item${active === 'start' ? ' on' : ''}`, href: '#/start' }, 'Start'),
    el('div', null,
      el('div', { class: 'rail-title' }, 'The route'),
      el('div', { class: 'route' }, ROUTE.map(([id, title, hint]) => el('a', {
        class: `stop${active === id ? ' on' : ''}`, href: `#/${id}`,
      }, el('i'), el('span', null, el('b', null, title), el('small', null, hint)))))),
    el('div', null,
      el('div', { class: 'rail-title' }, 'The bench'),
      el('a', { class: `rail-item${active === 'references' ? ' on' : ''}`, href: '#/references' }, 'References'),
      el('a', { class: `rail-item${active === 'tools' ? ' on' : ''}`, href: '#/tools' }, 'Tools',
        found === null ? null : el('span', { class: 'count' }, `${found} of ${tools.length}`))),
    el('div', { class: 'local' }, 'Stays on this computer', el('small', null, location.host)));
}

/* The plastome strip: the four regions to scale, with the genes on their strands. ---------------- */

function strip(gate, { labels = true, bracket = true } = {}) {
  const regions = gate.regions || [];
  const total = regions.reduce((sum, region) => sum + region.length, 0);
  const box = el('div', { class: 'strip' });
  if (!total) return box;
  const drawing = svg('svg', { viewBox: '0 0 1000 44', preserveAspectRatio: 'none', role: 'img',
    'aria-label': regions.map((r) => `${r.name} ${r.length} bp`).join(', ') });
  let left = 0;
  for (const region of regions) {
    const width = region.length * 1000 / total;
    drawing.append(svg('rect', { x: left, y: 15, width, height: 14, class: `r-${region.name}` }));
    left += width;
  }
  for (const [start, width, strand] of gate.ticks || []) {
    drawing.append(svg('rect', { class: 't', x: start * 10, width: Math.max(width * 10, 1), y: strand ? 2 : 31, height: 11 }));
  }
  if (bracket && gate.status && regions.length === 4) {
    const first = (regions[0].length + regions[1].length / 2) * 100 / total;
    const last = 100 - regions[3].length / 2 * 100 / total;
    const line = el('div');
    const label = el('span', null, el('b', null,
      `${gate.identity_pct.toFixed(2)} % identical, ${gate.status === 'passed' ? 'passed' : 'refused'}`));
    for (const piece of [line, label]) { piece.style.left = `${first}%`; piece.style.width = `${last - first}%`; }
    box.append(el('div', { class: `bracket ${gate.status}` }, line, label));
  }
  box.append(drawing);
  if (labels) {
    const row = el('div', { class: 'strip-labels' });
    for (const region of regions) {
      const cell = el('span', null, el('b', null, region.name), ' ', el('code', null, number(region.length)));
      cell.style.width = `${region.length * 100 / total}%`;
      row.append(cell);
    }
    box.append(row);
  }
  return box;
}

/* Choosing a file from the workspace. The server reads the disk, so a path survives Docker. ------- */

function chooseFile({ title = 'Choose a file', start = '', folders = false } = {}) {
  return new Promise((resolve) => {
    const list = el('div', { class: 'picker-list' });
    const where = el('div', { class: 'small mono' });
    const dialog = el('dialog', null,
      el('div', { class: 'card' }, el('h2', null, title), where, list,
        el('div', { class: 'row' },
          folders ? el('button', { class: 'btn go', type: 'button', onclick: () => done(current) }, 'Use this folder') : null,
          el('button', { class: 'btn', type: 'button', onclick: () => done(null) }, 'Close'))));
    let current = start;
    const done = (value) => { dialog.close(); dialog.remove(); resolve(value); };
    const open = async (path) => {
      const listing = await get('/api/browse', { path });
      current = listing.path;
      where.textContent = `/${listing.path}`;
      const up = listing.path.includes('/') ? listing.path.replace(/\/[^/]*$/, '') : '';
      list.replaceChildren(
        listing.path ? el('button', { class: 'folder', type: 'button', onclick: () => open(up) }, '.. up one folder') : '',
        ...listing.entries.map((entry) => el('button', {
          type: 'button', class: entry.folder ? 'folder' : '',
          onclick: () => (entry.folder ? open(entry.path) : done(entry.path)),
        }, entry.folder ? `${entry.name}/` : entry.name)));
    };
    document.body.append(dialog);
    dialog.showModal();
    open(start).catch(() => open(''));
  });
}

function pathField(labelText, value, change, options) {
  const input = el('input', { type: 'text', value: value || '', oninput: (event) => change(event.target.value) });
  return el('label', { class: 'field' }, labelText,
    el('span', { class: 'path-field' }, input,
      el('button', { class: 'btn', type: 'button', onclick: async (event) => {
        event.preventDefault();
        const chosen = await chooseFile({ title: labelText, ...options });
        if (chosen !== null) { input.value = chosen; change(chosen); }
      } }, 'Choose')));
}

function textField(labelText, value, change, note) {
  return el('label', { class: 'field' }, labelText,
    el('input', { type: 'text', value: value || '', oninput: (event) => change(event.target.value) }),
    note ? el('span', { class: 'tiny' }, note) : null);
}

function checkField(labelText, value, change) {
  return el('label', { class: 'check' }, el('input', { type: 'checkbox', checked: value, onchange: (event) => change(event.target.checked) }), labelText);
}

function selectField(labelText, options, value, change) {
  return el('label', { class: 'field' }, labelText, el('select', { onchange: (event) => change(event.target.value) },
    options.map((option) => el('option', { value: option, selected: option === value }, option))));
}

// Comma-separated words in a box, a list in the form.
function listField(labelText, change, note) {
  return textField(labelText, '', (value) => change(value.split(',').map((item) => item.trim()).filter(Boolean)), note);
}

// What most people never change, closed until someone wants it.
function more(title, ...fields) {
  return el('details', { class: 'more' }, el('summary', null, title), el('div', { class: 'field-grid' }, ...fields));
}

// A button that starts a command. What the server refuses is said beside the button.
function jobButton(command, build, words, kind = 'btn') {
  const button = el('button', { class: kind, type: 'button', onclick: async () => {
    if (button.nextElementSibling && button.nextElementSibling.matches('.note-warn')) button.nextElementSibling.remove();
    try {
      const job = await post('/api/jobs', { command, form: build() });
      location.hash = `#/job/${job.id}`;
    } catch (error) { button.after(fail(error)); }
  } }, words);
  return button;
}

/* The line a form will run, shown beside it, and the button that runs exactly that line. */

function commandBox(command, form) {
  const line = el('pre', { class: 'command' }, '');
  const start = el('button', { class: 'btn go', type: 'button', disabled: true, onclick: async () => {
    start.disabled = true;
    try {
      const job = await post('/api/jobs', { command, form: form() });
      location.hash = `#/job/${job.id}`;
    } catch (error) { line.textContent = String(error.message); }
  } }, 'Start');
  let asked = 0;
  const refresh = async (blocked) => {
    asked += 1;
    const mine = asked;
    const answer = await post('/api/command', { command, form: form() });
    if (mine !== asked) return;  // a newer keystroke has its own answer coming
    line.textContent = answer.line || answer.incomplete;
    start.disabled = Boolean(blocked) || !answer.line;
  };
  return { line, start, refresh };
}

/* Start ------------------------------------------------------------------------------------------ */

function doorArt(kind) {
  const art = svg('svg', { viewBox: '0 0 220 58', 'aria-hidden': 'true' });
  if (kind === 'reads') {
    for (const [x, y, w] of [[0, 8, 34], [44, 8, 26], [96, 8, 38], [156, 8, 30], [14, 20, 40], [70, 20, 34], [118, 20, 26], [176, 20, 40],
      [4, 32, 28], [46, 32, 38], [104, 32, 36], [150, 32, 24], [24, 44, 36], [78, 44, 28], [128, 44, 40], [182, 44, 30]]) {
      art.append(svg('rect', { x, y, width: w, height: 4, fill: '#13231D' }));
    }
  } else if (kind === 'sequence') {
    art.append(svg('rect', { x: 0, y: 23, width: 220, height: 12, rx: 2, fill: '#13231D' }));
  } else {
    const rows = kind === 'several' ? [10, 34] : [22];
    for (const y of rows) {
      [[0, 119, 'r-LSC'], [119, 38, 'r-IRa'], [157, 26, 'r-SSC'], [183, 37, 'r-IRb']].forEach(([x, width, name]) => {
        art.append(svg('rect', { x, y, width, height: 12, class: name }));
      });
    }
  }
  return art;
}

async function startPage() {
  const doors = [
    ['assemble', 'reads', 'Reads', 'Paired short reads or HiFi reads, FASTQ. Assemble a plastome and check it against the reads.', 'Assemble →'],
    ['annotate', 'sequence', 'A sequence', 'FASTA with no genes marked. Annotate it from reference records you choose.', 'Annotate →'],
    ['records', 'record', 'An annotated record', 'GenBank, one genome. Pass the IR gate, then choose the analyses.', 'Check the record →'],
    ['records', 'several', 'Several records', 'GenBank, two genomes or more. Compare diversity, Ka/Ks, synteny and phylogeny.', 'Compare →'],
  ];
  const latest = el('tbody');
  const machine = el('div', { class: 'card' }, 'Asking the tools what they are…');
  page('start',
    el('header', null, el('h1', null, 'What do you have?'),
      el('p', { class: 'lead' }, 'Every route ends in the same tables, figures and report. Your data stays on this computer.')),
    el('section', { class: 'grid', 'aria-label': 'Starting points' }, doors.map(([route, art, title, text, cta]) => el('a', { class: 'door', href: `#/${route}` },
      doorArt(art), el('b', null, title), el('p', null, text), el('em', null, cta)))),
    el('div', { class: 'columns' },
      el('section', { class: 'grow' },
        el('div', { class: 'row between' }, el('h2', null, 'Latest runs'), el('a', { href: '#/runs' }, 'All runs')),
        el('div', { class: 'table-box' }, el('table', null,
          el('thead', null, el('tr', null, ['Status', 'Dataset', 'Analyses', 'Started'].map((name) => el('th', { scope: 'col' }, name)))), latest))),
      el('section', { class: 'side' }, el('h2', null, 'This computer'), machine)));
  const runs = (await get('/api/runs')).runs.slice(0, 5);
  latest.replaceChildren(...(runs.length ? runs.map(runRow) : [el('tr', null, el('td', { colspan: 4, class: 'empty' }, 'No run yet.'))]));
  const tools = await loadTools();
  const missing = tools.tools.filter((tool) => tool.status === 'missing').map((tool) => tool.name);
  machine.replaceChildren(
    el('div', { class: 'row between' }, el('span', { class: 'small' }, 'Your folder'), el('code', { class: 'small' }, session.about.workspace)),
    el('div', { class: 'row between' }, el('span', { class: 'small' }, 'Tools found'),
      el('a', { href: '#/tools' }, `${tools.tools.length - missing.length} of ${tools.tools.length}`)),
    missing.length ? el('p', { class: 'small' }, `Not found: ${missing.join(', ')}. Analyses that need one say so before they start.`) : null);
}

async function loadTools() {
  if (!session.tools) session.tools = await get('/api/tools');
  return session.tools;
}

function runRow(run) {
  return el('tr', null,
    el('td', null, word(run.status)),
    el('td', null, el('a', { href: `#/run/${run.folder}/${run.run_id}` }, run.dataset)),
    el('td', { class: 'num' }, run.total ? `${run.done} of ${run.total}${run.kind === 'batch' ? ' samples' : ''}` : 'figures'),
    el('td', { class: 'num' }, when(run.started_at)));
}

/* Records ---------------------------------------------------------------------------------------- */

function gateOf(record) {
  const gates = record.gates || {};
  return session.viaSequence.has(record.path) ? gates.sequence : (gates.annotation || gates.sequence);
}

function recordCard(record, refreshTop) {
  if (record.unreadable) {
    return el('article', { class: 'card' }, el('h2', null, record.path), el('div', { class: 'note-warn' }, record.unreadable));
  }
  const gate = gateOf(record) || {};
  const genes = record.genes.total;
  const facts = [`${number(record.length)} bp`, `GC ${record.gc_pct.toFixed(2)} %`, genes ? `${genes} genes` : 'no genes marked', record.path.split('/').pop()];
  const refused = gate.status === 'refused' && !session.viaSequence.has(record.path);
  const card = el('article', { class: 'card' },
    el('div', { class: 'row between' },
      el('label', { class: 'check' },
        record.annotated ? el('input', { type: 'checkbox', checked: session.chosen.has(record.path), onchange: (event) => {
          if (event.target.checked) session.chosen.add(record.path); else session.chosen.delete(record.path);
          refreshTop();
        } }) : null,
        el('span', null, el('b', null, species(record.organism), ' ', el('code', { class: 'small' }, record.accession || '')),
          el('span', { class: 'small', style: null }, el('br'), facts.join(' · ')))),
      record.annotated
        ? el('span', { class: `pill ${gate.status === 'passed' ? 'ok' : 'warn'}` },
          gate.status === 'passed' ? (session.viaSequence.has(record.path) ? 'Gate passed, IR read from the sequence' : 'Gate passed') : 'Refused at the gate')
        : el('a', { class: 'btn', href: `#/annotate?input=${encodeURIComponent(record.path)}` }, 'Annotate first')),
    strip(gate, { bracket: record.annotated }));
  if (refused) {
    const annotated = gate.annotated_lengths.map(number).join(' and ');
    const other = record.gates.sequence || {};
    card.append(el('div', { class: 'note-warn' },
      el('p', null, `The record marks its IR copies as ${annotated} bp, which are ${gate.identity_pct.toFixed(2)} % identical. `,
        gate.denovo_length ? `The sequence itself holds a repeat of ${number(gate.denovo_length)} bp. ` : '',
        'Position-based analyses cannot run on the marked pair.'),
      el('div', { class: 'row' },
        other.status === 'passed' ? el('button', { class: 'btn way', type: 'button', onclick: () => {
          session.viaSequence.add(record.path);
          card.replaceWith(recordCard(record, refreshTop));
        } }, 'Read the IR from the sequence') : null,
        el('span', { class: 'small' }, 'Or correct the record’s repeat_region marks and add it again. The gate has no off switch.'))));
  }
  if (session.viaSequence.has(record.path)) {
    card.append(el('div', { class: 'row between note-ok' },
      el('span', null, `${number(gate.length)} bp, ${gate.identity_pct.toFixed(2)} % identical. Every run of this record will say the IR came from the sequence.`),
      el('button', { class: 'link', type: 'button', onclick: () => {
        session.viaSequence.delete(record.path);
        card.replaceWith(recordCard(record, refreshTop));
      } }, 'Go back to the marked pair')));
  }
  return card;
}

async function recordsPage() {
  const list = el('div', { class: 'grow' }, 'Reading the records…');
  const next = el('a', { class: 'btn go', href: '#/analyse' }, 'Analyse');
  const refreshTop = () => {
    const count = session.chosen.size;
    next.textContent = count > 1 ? `Compare these ${count}` : (count === 1 ? 'Analyse this record' : 'Analyse');
  };
  const upload = el('input', { type: 'file', multiple: true, hidden: true, onchange: async (event) => {
    for (const file of event.target.files) {
      await call(`/api/upload?${new URLSearchParams({ name: file.name, to: 'records' })}`, {
        method: 'POST', headers: { 'Content-Type': 'application/octet-stream' }, body: file });
    }
    recordsPage();
  } });
  page('records',
    el('header', { class: 'page-head' },
      el('div', null, el('h1', null, 'Records'),
        el('p', { class: 'lead' }, 'A record passes the IR gate before any analysis that reads positions. The two inverted-repeat copies must be at least 99.9 % identical.')),
      el('div', { class: 'row' }, upload,
        el('button', { class: 'btn', type: 'button', onclick: () => upload.click() }, 'Add records'), next)),
    list);
  refreshTop();
  const { records } = await get('/api/records');
  if (!records.length) {
    list.replaceChildren(el('div', { class: 'card empty' }, 'No record yet. Add a GenBank or FASTA file, or put it in the records folder.'));
    return;
  }
  const places = records.map((entry) => el('article', { class: 'card' }, el('h2', null, entry.name), el('p', { class: 'small' }, 'Reading the record and asking the IR gate…')));
  list.replaceChildren(...places);
  await Promise.all(records.map(async (entry, index) => {
    places[index].replaceWith(recordCard(await get('/api/record', { path: entry.path }), refreshTop));
  }));
}

/* A parameter's field, built from what the module declares. ------------------------------------- */

function paramField(param, values, changed) {
  const current = values[param.name] ?? param.default ?? '';
  const change = (value) => {
    if (String(value) === String(param.default ?? '')) delete values[param.name]; else values[param.name] = value;
    changed();
  };
  let control;
  if (param.choices && param.choices.length) {
    control = el('select', { onchange: (event) => change(event.target.value) },
      param.choices.map((choice) => el('option', { value: choice, selected: String(choice) === String(current) }, choice)));
  } else if (param.type === 'bool') {
    control = el('label', { class: 'check' },
      el('input', { type: 'checkbox', checked: current === true || current === 'true' || current === 'yes',
        onchange: (event) => change(event.target.checked ? 'yes' : 'no') }), 'yes');
  } else if (param.type === 'file' || param.type === 'dir') {
    return pathField(param.name, String(current), change, { folders: param.type === 'dir' });
  } else {
    control = el('input', { type: 'text', value: String(current), oninput: (event) => {
      event.target.classList.toggle('changed', event.target.value !== String(param.default ?? ''));
      change(event.target.value);
    } });
  }
  return el('div', { class: 'field' },
    el('span', { class: 'param-name' }, param.name), control,
    el('span', { class: 'small' }, param.help),
    el('span', { class: 'tiny' }, param.default === null || param.default === undefined ? 'Default: not set' : `Default: ${param.default}`));
}

/* Analyse ---------------------------------------------------------------------------------------- */

async function analysePage() {
  if (!session.modules) session.modules = await get('/api/modules');
  const tools = await loadTools();
  const missing = new Set(tools.tools.filter((tool) => tool.status === 'missing').map((tool) => tool.name));
  const all = session.modules.modules.filter((module) => module.group !== 'mito' && !module.command.includes('annotate'));
  const inputs = [...session.chosen];
  const form = session.form = session.form || { mode: inputs.length > 1 ? 'compare' : 'one', on: new Set(), focus: null, sets: {}, label: '', dataset: '', figures: tools.extras.figures, docx: false, no_cache: false };
  if (inputs.length > 1 && form.mode === 'one') { form.mode = 'compare'; form.on = new Set(); }
  const place = (module) => { const at = ORDER.indexOf(module.id); return at < 0 ? ORDER.length : at; };
  const shown = () => all
    .filter((module) => (form.mode === 'compare' ? module.genomes.min > 1 : module.in_batch && module.genomes.min === 1))
    .sort((one, other) => place(one) - place(other));
  const records = await Promise.all(inputs.map((path) => get('/api/record', { path })));
  const needsTool = (module) => module.tools.find((tool) => missing.has(tool));
  const viaSequence = () => inputs.some((path) => session.viaSequence.has(path));
  const blockedByGate = () => records.some((record) => {
    const gate = gateOf(record);
    return gate && gate.status !== 'passed';
  });

  const list = el('section', { class: 'card pick', 'aria-label': 'Analyses' });
  const parameters = el('section', { class: 'card', 'aria-label': 'Parameters' });
  const plan = el('ol', { class: 'plan' });
  const notice = el('div');
  const payload = () => {
    const sets = { ...form.sets };
    if (viaSequence()) sets.ir_source = 'sequence';
    return { modules: shown().filter((module) => form.on.has(module.id)).map((module) => module.id), inputs, one_by_one: form.mode === 'batch', label: form.label, dataset: form.dataset, groups: form.groups, figures: form.figures, docx: form.docx, no_cache: form.no_cache, sets };
  };
  const box = commandBox('run', payload);
  const redraw = () => {
    const modules = shown();
    list.replaceChildren(
      el('div', { class: 'row between small' }, el('b', null, form.mode === 'compare' ? 'Two genomes or more' : 'One genome each'), el('span', { class: 'mono' }, `${[...form.on].length} chosen`)),
      ...modules.map((module) => {
        const lacking = needsTool(module);
        return el('div', { class: `choice${form.focus === module.id ? ' on' : ''}` },
          el('input', { type: 'checkbox', 'aria-label': module.title, checked: form.on.has(module.id), onchange: (event) => {
            if (event.target.checked) form.on.add(module.id); else form.on.delete(module.id);
            form.focus = module.id; redraw();
          } }),
          el('button', { type: 'button', onclick: () => { form.focus = module.id; redraw(); } }, module.title,
            lacking ? el('small', { class: 'warn' }, `needs ${lacking}, not found`) : (module.ir_gate ? el('small', null, 'passes the IR gate first') : null)));
      }));
    const focus = modules.find((module) => module.id === form.focus) || modules[0];
    form.focus = focus ? focus.id : null;
    parameters.replaceChildren(...(focus ? [
      el('h2', null, focus.title),
      el('p', { class: 'small' }, [focus.genomes.min > 1 ? 'Two genomes or more' : 'One genome', focus.ir_gate ? 'Passes the IR gate first' : 'Does not read IR positions', focus.tools.length ? `Uses ${focus.tools.join(', ')}` : 'No external tool'].join(' · ')),
      focus.reference_param ? pathField('The record you chose as the naming standard', form.sets[focus.reference_param], (value) => {
        if (value) form.sets[focus.reference_param] = value; else delete form.sets[focus.reference_param];
        box.refresh(blockers());
      }, { start: 'references' }) : null,
      focus.params.length ? el('div', { class: 'field-grid' }, focus.params.filter((param) => !param.flag).map((param) => paramField(param, form.sets, () => box.refresh(blockers()))))
        : el('p', { class: 'small' }, 'This analysis has no parameters.'),
    ] : []));
    const chosen = modules.filter((module) => form.on.has(module.id));
    const steps = [];
    if (chosen.some((module) => module.ir_gate)) {
      const gate = records.map(gateOf).find((item) => item && item.status !== 'passed');
      steps.push(gate ? ['Refused', 'warn', 'IR gate', `${gate.identity_pct} % identical, the gate needs 99.90 %`]
        : ['Passed', 'ok', 'IR gate', viaSequence() ? 'read from the sequence' : 'read from the annotation']);
    }
    for (const module of chosen) {
      const lacking = needsTool(module);
      steps.push(lacking ? ['Needs tool', 'warn', module.title, `${lacking} was not found`] : ['Ready', 'ok', module.title, module.tools.join(', ') || 'no external tool']);
    }
    plan.replaceChildren(...steps.map(([label, tone, title, detail]) => el('li', null, el('b', { class: `word-${tone}` }, label), el('span', null, title, el('small', null, detail)))));
    const reason = blockers();
    notice.replaceChildren(...(reason ? [el('div', { class: 'note-warn' }, reason.words, reason.link ? el('a', { href: reason.link }, reason.linkWords) : null)] : []));
    box.refresh(reason);
  };
  const blockers = () => {
    const chosen = shown().filter((module) => form.on.has(module.id));
    if (!inputs.length) return { words: 'Choose a record first.', link: '#/records', linkWords: 'Open Records' };
    if (!chosen.length) return { words: 'Choose at least one analysis.' };
    if (chosen.some((module) => module.ir_gate) && blockedByGate()) return { words: 'A record stopped at the IR gate. Read its IR from the sequence, or leave out the analyses that pass the gate.', link: '#/records', linkWords: 'Open Records' };
    const lacking = chosen.find(needsTool);
    if (lacking) return { words: `${lacking.title} needs ${needsTool(lacking)}, which this computer does not have.`, link: '#/tools', linkWords: 'Open Tools' };
    return null;
  };
  const setMode = (mode) => { form.mode = mode; form.on = new Set(); form.focus = null; analysePage(); };
  const mode = ([id, words]) => el('button', { type: 'button', 'aria-pressed': String(form.mode === id), onclick: () => setMode(id) }, words);
  page('analyse',
    el('header', { class: 'page-head' },
      el('div', null, el('h1', null, 'Analyse'), el('p', { class: 'lead' }, 'Choose the analyses. The plan on the right shows what will run before anything starts.')),
      el('div', { class: 'seg', role: 'group', 'aria-label': 'How many genomes' }, [['one', 'One genome'], ['compare', 'Compare genomes'], ['batch', 'Many, one by one']].map(mode))),
    el('div', { class: 'columns' }, list,
      el('div', { class: 'grow' },
        el('section', { class: 'card', 'aria-label': 'Run details' },
          el('div', null, el('div', { class: 'small' }, el('b', null, 'Records')),
            el('div', { class: 'row' }, records.map((record) => el('span', { class: 'note-plain' }, species(record.organism), ' ', el('code', { class: 'small' }, record.accession))), el('a', { href: '#/records' }, inputs.length ? 'Change' : 'Choose records'))),
          el('div', { class: 'field-grid' },
            textField('Where the data came from', form.label, (value) => { form.label = value; box.refresh(blockers()); }, 'Printed with every table.'),
            textField('Dataset name', form.dataset, (value) => { form.dataset = value; box.refresh(blockers()); }, 'The folder inside runs.'),
            form.mode === 'batch' ? null : pathField('Sample groups (optional)', form.groups, (value) => { form.groups = value; box.refresh(blockers()); }, { start: 'records' })),
          el('div', { class: 'row' }, [['figures', 'Draw figures'], ['docx', 'Word reports']].map(([key, words]) => el('label', { class: 'check' }, el('input', { type: 'checkbox', checked: form[key], disabled: !tools.extras[key], onchange: (event) => { form[key] = event.target.checked; box.refresh(blockers()); } }), tools.extras[key] ? words : `${words}: needs plastid-suite[${key}]`)),
            el('label', { class: 'check' }, el('input', { type: 'checkbox', checked: !form.no_cache, onchange: (event) => { form.no_cache = !event.target.checked; box.refresh(blockers()); } }), 'Reuse stored results'))),
        parameters),
      el('aside', { class: 'card side', 'aria-label': 'Plan' }, el('h2', null, 'The plan'), plan, box.start, notice,
        el('div', { class: 'small' }, el('b', null, 'The same run on the command line')), box.line)));
  box.start.textContent = 'Run';
  redraw();
}

/* A command in progress -------------------------------------------------------------------------- */

async function jobPage(id) {
  const summary = el('section', { class: 'card summary' });
  const steps = el('section', { class: 'card pick' });
  const log = el('pre', { class: 'log grow' });
  const actions = el('div', { class: 'row' });
  const head = el('h1', null, 'Working');
  const note = el('p', { class: 'lead' }, 'Closing this tab does not stop the work. Come back at any time.');
  page('runs',
    el('header', { class: 'page-head' }, el('div', null, head, note), actions), summary,
    el('div', { class: 'columns' }, steps, log));
  let stopped = false;
  const refresh = async () => {
    if (stopped) return;
    const job = await get(`/api/jobs/${id}`);
    const run = job.run ? await get(`/api/run/${job.run.dataset}/${job.run.run_id}`) : null;
    log.textContent = job.output || 'Starting…';
    summary.replaceChildren(
      el('div', null, el('span', null, 'Status'), el('b', null, word(job.status))),
      run && run.total ? el('div', null, el('span', null, run.kind === 'batch' ? 'Samples' : 'Analyses'), el('b', null, `${run.done} of ${run.total}`)) : null,
      el('div', { class: 'wide' }, el('span', null, 'Command'), el('b', null, job.line)));
    if (run) {
      const parts = run.kind === 'batch' ? run.samples.map((sample) => ({ id: sample.sample_id, status: sample.status })) : run.steps;
      const rows = parts.map((step) => el('div', { class: `step-row ${['ok', 'cached'].includes(step.status) ? 'done' : (step.status === 'failed' ? 'failed' : '')}` }, el('i'), step.id, el('small', null, step.status || '')));
      if (job.running) rows.push(el('div', { class: 'step-row now' }, el('i'), run.current || 'starting', el('small', null, 'running')));
      steps.replaceChildren(...rows);
    }
    actions.replaceChildren(...(job.running && job.run ? [el('button', { class: 'btn stop-run', type: 'button', onclick: async () => {
      await post(`/api/run/${job.run.dataset}/${job.run.run_id}/cancel`);
    } }, 'Cancel run')] : []));
    if (!job.running) {
      stopped = true;
      for (const timer of timers) clearInterval(timer);
      head.textContent = statusOf(job.status)[0];
      note.textContent = job.status === 'cancelled' ? 'What had finished stays in the run directory.' : '';
      actions.replaceChildren(job.run
        ? el('a', { class: 'btn go', href: `#/run/${job.run.dataset}/${job.run.run_id}` }, 'Open the results')
        : el('button', { class: 'btn', type: 'button', onclick: () => history.back() }, 'Back'));
    }
  };
  every(1500, refresh);
  refresh();
}

/* Results ---------------------------------------------------------------------------------------- */

async function runsPage() {
  const { runs } = await get('/api/runs');
  const filter = { show: 'all' };
  const body = el('tbody');
  const buttons = el('div', { class: 'row' });
  const redraw = () => {
    const statuses = ['all', ...new Set(runs.map((run) => run.status))];
    buttons.replaceChildren(...statuses.map((status) => el('button', { class: 'btn', type: 'button', 'aria-pressed': String(filter.show === status), onclick: () => { filter.show = status; redraw(); } },
      `${status === 'all' ? 'All' : statusOf(status)[0]} ${status === 'all' ? runs.length : runs.filter((run) => run.status === status).length}`)));
    const rows = runs.filter((run) => filter.show === 'all' || run.status === filter.show);
    body.replaceChildren(...(rows.length ? rows.map(runRow) : [el('tr', null, el('td', { colspan: 4, class: 'empty' }, 'No run yet.'))]));
  };
  page('runs',
    el('header', null, el('h1', null, 'Results'), el('p', { class: 'lead' }, 'Every run in this workspace, newest first. Runs started on the command line appear here too.')),
    buttons,
    el('div', { class: 'table-box' }, el('table', null, el('thead', null, el('tr', null, ['Status', 'Dataset', 'Analyses', 'Started'].map((name) => el('th', { scope: 'col' }, name)))), body)));
  redraw();
  every(5000, async () => {
    const fresh = (await get('/api/runs')).runs;
    if (fresh.length !== runs.length || fresh.some((run, index) => run.status !== runs[index].status)) runsPage();
  });
}

function provenance(run) {
  const bits = [];
  if (run.inputs.length === 1) bits.push(`Read from ${run.inputs[0].path.split('/').pop()}${run.inputs[0].sha256 ? ` (sha256 ${run.inputs[0].sha256.slice(0, 8)}…)` : ''}.`);
  else if (run.inputs.length) bits.push(`Read from ${run.inputs.length} records, listed under How this was made.`);
  if (run.ir && run.ir.source) bits.push(`IR read from ${run.ir.annotation_used ? 'the annotation' : 'the sequence'}.`);
  const label = run.params.data_label || run.data_label;
  if (label) bits.push(`Data from ${label}.`);
  else bits.push('Data source: not labelled.');
  if (run.package_version) bits.push(`PlastidToolSuite ${run.package_version}.`);
  return bits.join(' ');
}

async function tableView(entry) {
  const table = await get('/api/table', { path: entry.path });
  const notice = table.total > table.rows.length ? el('p', { class: 'tiny' }, `Showing ${table.rows.length} of ${number(table.total)} rows. The file holds them all.`) : null;
  return el('div', { class: 'card' },
    el('div', { class: 'row between' }, el('h3', null, table.title || entry.name), el('a', { href: `/api/file?${new URLSearchParams({ path: entry.path, download: 1 })}` }, 'CSV')),
    el('div', { class: 'table-box' }, el('table', null,
      el('thead', null, el('tr', null, table.columns.map((column) => el('th', { scope: 'col' }, column)))),
      el('tbody', null, table.rows.map((row) => el('tr', null, row.map((cell) => (cell === ''
        ? el('td', { class: 'notmeasured' }, 'not measured')
        : el('td', { class: /^-?[0-9][0-9,.]*$/.test(cell) ? 'num' : '' }, cell)))))))),
    notice,
    el('p', { class: 'small' }, el('b', null, 'Scope. '), table.scope || (table.columns.includes('scope') ? 'each row carries its own, in the scope column' : 'not stated'), table.coordinates ? ` | coordinates: ${table.coordinates}` : ''));
}

function figureCard(figure, label, path, pick, pressed) {
  return el('button', { class: 'figure-card', type: 'button', 'aria-pressed': String(pressed), onclick: pick },
    el('span', { class: 'frame' }, el('img', { src: `/api/preview?${new URLSearchParams({ path })}`, alt: '', loading: 'lazy' })),
    el('span', null, el('b', null, label.replaceAll('_', ' ')), el('br'), el('small', null, figure.id)));
}

async function figuresView(run) {
  const figures = run.figures;
  const views = figures.drawn.flatMap((figure) => Object.entries(figure.views).map(([label, path]) => ({ figure, label, path })));
  const big = el('section', { class: 'card' });
  const grid = el('section', { class: 'grid gallery', 'aria-label': 'Gallery' });
  let open = views[0];
  const redraw = async () => {
    grid.replaceChildren(...views.map((view) => figureCard(view.figure, view.label, view.path, () => { open = view; redraw(); }, view === open)));
    if (!open) { big.replaceChildren(el('p', { class: 'empty' }, 'This run holds no figure yet.')); return; }
    const caption = open.figure.caption ? (await get('/api/text', { path: open.figure.caption })).text : '';
    big.replaceChildren(el('div', { class: 'figure-big' },
      el('div', null, el('img', { src: `/api/preview?${new URLSearchParams({ path: open.path })}`, alt: open.label })),
      el('div', null, el('h2', null, open.label.replaceAll('_', ' ')),
        el('div', { class: 'row' }, open.figure.files.map((file) => el('a', { class: 'btn', href: `/api/file?${new URLSearchParams({ path: file.path, download: 1 })}` }, file.name.split('.').pop().toUpperCase()))),
        caption ? el('pre', { class: 'small' }, caption) : null)));
  };
  const absent = [
    figures.error ? el('li', null, el('b', { class: 'word-bad' }, 'Failed'), el('span', null, 'All figures', el('small', null, figures.error))) : null,
    ...figures.failed.map((entry) => el('li', null, el('b', { class: 'word-bad' }, 'Failed'), el('span', null, entry.figure, el('small', null, entry.reason)))),
    ...figures.skipped.map((entry) => el('li', null, el('b', { class: 'word-quiet' }, 'Not drawn'), el('span', null, entry.figure, el('small', null, entry.reason)))),
    ...figures.missing_renderers.map((entry) => el('li', null, el('b', { class: 'word-tool' }, 'Missing'), el('span', null, entry.renderer, el('small', null, entry.reason)))),
  ];
  await redraw();
  return [big, grid, absent.filter(Boolean).length ? el('section', { class: 'card' }, el('h2', null, 'What this run did not draw'), el('ol', { class: 'plan' }, absent)) : null];
}

function samplesView(run) {
  return el('section', { class: 'card' }, el('h2', null, 'Samples'),
    el('p', { class: 'small' }, 'One run each. A failed sample is listed with its reason and adds no number to the table for all.'),
    el('div', { class: 'table-box' }, el('table', null,
      el('thead', null, el('tr', null, ['Status', 'Species', 'Accession', 'IR gate', 'Its run'].map((name) => el('th', { scope: 'col' }, name)))),
      el('tbody', null, run.samples.map((sample) => [
        el('tr', null,
          el('td', null, word(sample.status)),
          el('td', null, species(sample.species)),
          el('td', { class: 'num' }, sample.accession),
          el('td', null, sample.ir_gate || 'not measured'),
          el('td', null, sample.run ? el('a', { href: `#/run/${sample.run.dataset}/${sample.run.run_id}` }, 'Open') : '')),
        sample.error ? el('tr', null, el('td', { colspan: 5, class: 'small' }, sample.error)) : null,
      ])))));
}

// A run's records on the same scale. A record that has left the workspace is simply not drawn.
async function drawMaps(run, box) {
  const answers = await Promise.all(run.inputs.slice(0, MAPS_SHOWN).map((input) => get('/api/record', { path: input.path }).catch(() => null)));
  const source = run.params.ir_source === 'sequence' ? 'sequence' : 'annotation';
  const drawn = answers.filter(Boolean).map((record) => {
    const gate = record.gates ? (record.gates[source] || record.gates.sequence) : null;
    if (!gate || !gate.regions) return null;
    return el('div', null, el('div', { class: 'small' }, species(record.organism), ` ${number(record.length)} bp`), strip(gate, { bracket: answers.length === 1 }));
  }).filter(Boolean);
  if (drawn.length) { box.replaceChildren(...drawn); box.hidden = false; }
}

async function runPage(dataset, runId, open) {
  const run = await get(`/api/run/${dataset}/${runId}`);
  if (!session.modules) session.modules = await get('/api/modules');
  const titled = (id) => (session.modules.modules.find((module) => module.id === id) || { title: id }).title;
  const tables = run.steps.flatMap((step) => step.tables.map((table) => ({ ...table, step: step.id })));
  const batch = run.kind === 'batch';
  if (!open && batch) open = 'samples';
  if (!open && !run.steps.length) open = 'figures';
  const body = el('section', { class: 'grow' });
  const nav = el('nav', { class: 'nav-list', 'aria-label': 'Sections of this run' });
  const go = (section) => { location.hash = `#/run/${dataset}/${runId}/${section}`; };
  const entry = (id, title) => el('button', { type: 'button', 'aria-current': String(open === id), onclick: () => go(id) }, title);
  nav.append(...[
    batch ? el('h3', null, 'This batch') : null, batch ? entry('samples', 'Samples') : null, run.aggregate ? entry('aggregate', 'One table for all') : null,
    run.steps.length ? el('h3', null, 'Analyses') : null, ...run.steps.map((step) => entry(step.id, titled(step.id))),
    el('h3', null, 'This run'), entry('figures', 'Figures'), entry('how', 'How this was made'), entry('files', 'Files'),
  ].filter(Boolean));
  const tiles = el('section', { class: 'grid tiles', 'aria-label': 'Key numbers' });
  const maps = el('section', { class: 'card', 'aria-label': 'Plastome maps', hidden: true });
  page('runs',
    el('nav', { class: 'crumb', 'aria-label': 'Breadcrumb' }, el('a', { href: '#/runs' }, 'Results'), ` / ${dataset} / `, el('code', null, runId)),
    el('header', { class: 'page-head' },
      el('div', null, el('h1', null, dataset), el('p', { class: 'lead' }, provenance(run))),
      el('div', { class: 'row' }, pill(run.status),
        el('a', { class: 'btn', href: `#/run/${dataset}/${runId}/figures` }, 'All figures'),
        el('a', { class: 'btn', href: `#/export?run=${dataset}/${runId}` }, 'Export'),
        run.status === 'running' ? el('button', { class: 'btn stop-run', type: 'button', onclick: () => post(`/api/run/${dataset}/${runId}/cancel`) }, 'Cancel run') : null)),
    maps, tiles,
    run.status === 'refused' || run.status === 'failed' ? el('div', { class: 'note-warn' }, run.steps.map((step) => step.error).filter(Boolean).join(' ')) : null,
    el('div', { class: 'columns' }, nav, body));
  if (run.status === 'running') {
    every(3000, async () => {
      const fresh = await get(`/api/run/${dataset}/${runId}`);
      if (fresh.status !== run.status || fresh.done !== run.done) runPage(dataset, runId, open);
    });
  }
  drawMaps(run, maps);
  if (batch) tiles.append(el('button', { class: 'tile', type: 'button', onclick: () => go('samples') }, el('span', null, 'Samples'), el('b', null, `${run.done} of ${run.total}`), el('small', null, 'finished')));
  for (const [file, column, label, noteWords] of HIGHLIGHTS) {
    const table = tables.find((candidate) => candidate.name === file);
    if (!table) continue;
    const data = await get('/api/table', { path: table.path });
    const index = data.columns.indexOf(column);
    if (index < 0 || !data.rows.length) continue;
    tiles.append(el('button', { class: 'tile', type: 'button', onclick: () => go(table.step) }, el('span', null, label), el('b', null, data.rows[0][index] || 'not measured'), el('small', null, noteWords)));
  }
  if (open === 'figures') { body.replaceChildren(...(await figuresView(run)).filter(Boolean)); return; }
  if (open === 'samples') { body.append(samplesView(run)); return; }
  if (open === 'aggregate') { body.append(await tableView({ path: run.aggregate, name: 'aggregate.csv' })); return; }
  if (open === 'how') {
    body.append(el('section', { class: 'card' }, el('h2', null, 'How this was made'),
      el('p', { class: 'small' }, 'Every tool is recorded with its version, its command line and its exit code.'),
      el('pre', { class: 'command' }, (run.command || []).join(' ')),
      el('ul', { class: 'small' }, run.inputs.map((input) => el('li', null, input.path, ' ', input.sha256 ? el('code', null, `sha256 ${input.sha256.slice(0, 16)}…`) : null))),
      el('div', { class: 'table-box' }, el('table', null,
        el('thead', null, el('tr', null, ['Step', 'Tool', 'Version', 'Calls', 'Exit'].map((name) => el('th', { scope: 'col' }, name)))),
        el('tbody', null, run.steps.flatMap((step) => [...new Set(step.tools.map((tool) => tool.name))].map((name) => {
          const calls = step.tools.filter((tool) => tool.name === name);
          return el('tr', null, el('td', null, step.id), el('td', null, name), el('td', { class: 'num' }, calls[0].version), el('td', { class: 'num' }, calls.length), el('td', { class: 'num' }, [...new Set(calls.map((tool) => tool.rc))].join(', ')));
        })))))));
    return;
  }
  if (open === 'files') {
    const listing = await get('/api/browse', { path: run.path });
    body.append(el('section', { class: 'card' }, el('h2', null, 'Files'),
      el('p', { class: 'small' }, 'A run directory is never overwritten. This is the same folder the command line writes.'),
      el('pre', { class: 'command' }, run.path), el('ul', null, listing.entries.map((item) => el('li', { class: 'mono small' }, item.folder ? `${item.name}/` : item.name)))));
    return;
  }
  const step = run.steps.find((item) => item.id === open) || run.steps[0];
  if (!step) { body.append(el('div', { class: 'card empty' }, 'This run finished no step.')); return; }
  body.append(el('h2', null, titled(step.id)));
  if (step.error) body.append(el('div', { class: 'note-warn' }, step.error));
  for (const table of step.tables) body.append(await tableView(table));
  if (step.report) {
    const report = await get('/api/text', { path: step.report });
    body.append(el('section', { class: 'card' }, el('h3', null, 'Draft report'), el('pre', { class: 'small' }, report.text)));
  }
}

/* Export ----------------------------------------------------------------------------------------- */

async function exportPage(query) {
  const { runs } = await get('/api/runs');
  const finished = runs.filter((run) => !['running', 'refused'].includes(run.status));
  const address = (run) => `${run.folder}/${run.run_id}`;
  const initial = query.get('run') || (finished[0] ? address(finished[0]) : '');
  const current = finished.find((run) => address(run) === initial);
  const runPath = `${session.about.runs}/${initial}`;
  const action = (title, text, command, words) => el('div', { class: 'card' }, el('h2', null, title), el('p', { class: 'small' }, text),
    current && current.kind === 'batch'
      ? el('p', { class: 'small' }, 'This works on one run: open a sample of this batch, then export it.')
      : jobButton(command, () => ({ run: runPath }), words));
  const drawing = { figures: [], sets: {} };
  const drawBox = commandBox('figures', () => ({ run: runPath, figures: drawing.figures, sets: drawing.sets }));
  const refreshDraw = () => drawBox.refresh(!current);
  const catalogue = current ? await get('/api/figures', { run: runPath }) : { figures: [], options: [] };
  const junction = { inputs: [], label: '', dataset: 'ir-junction', cpjsdraw: '', tag: '', width: '', tree: '', order_from: '', drop_pseudo: false };
  const junctionBox = commandBox('ir-junction', () => junction);
  const refreshJunction = () => junctionBox.refresh();
  const records = (await get('/api/records')).records;
  page('export',
    el('header', { class: 'page-head' },
      el('div', null, el('h1', null, 'Export'), el('p', { class: 'lead' }, 'Build what leaves the toolkit from a finished run. Nothing is computed again.')),
      el('label', { class: 'field' }, 'Run', el('select', { onchange: (event) => { location.hash = `#/export?run=${event.target.value}`; } },
        finished.map((run) => el('option', { value: address(run), selected: address(run) === initial }, `${run.dataset}, ${when(run.started_at)}`))))),
    el('section', { class: 'grid' },
      action('Supplementary package', 'Tables, figures, reports and an index, in one folder.', 'package', 'Build the package'),
      action('Technical QC report', 'What ran, tool versions, inputs, audits and figures.', 'qc', 'Write the report'),
      el('div', { class: 'card' }, el('h2', null, 'NCBI submission'), el('p', { class: 'small' }, 'Write .fsa and .tbl from a checked record. With your .sbt template, table2asn adds a .sqn.'), el('a', { class: 'btn', href: '#/analyse' }, 'Set it up in Analyse'))),
    el('section', { class: 'card' }, el('h2', null, 'Figures'),
      el('p', { class: 'small' }, 'Leave every box empty to draw all the run’s tables can feed, or tick the ones to draw again. A publication figure comes as PNG and PDF, a quick one as SVG, EPS and PDF.'),
      el('ol', { class: 'plan' }, catalogue.figures.filter((figure) => figure.can_draw).map((figure) => el('li', null,
        el('input', { type: 'checkbox', 'aria-label': figure.id, onchange: (event) => {
          drawing.figures = event.target.checked ? [...drawing.figures, figure.id] : drawing.figures.filter((id) => id !== figure.id);
          refreshDraw();
        } }),
        el('span', null, figure.id, el('small', null, figure.drawn ? `drawn, ${figure.renderer === 'R' ? 'publication figure' : 'quick figure'}` : figure.reason))))),
      more('Drawing options', catalogue.options.map((param) => paramField(param, drawing.sets, refreshDraw))),
      el('details', { class: 'more' }, el('summary', null, 'What this run has no data for'),
        el('ol', { class: 'plan' }, catalogue.figures.filter((figure) => !figure.can_draw).map((figure) => el('li', null,
          el('b', { class: 'word-quiet' }, 'No data'), el('span', null, figure.id, el('small', null, figure.reason)))))),
      drawBox.line, drawBox.start),
    el('section', { class: 'card' }, el('h2', null, 'IR junctions across species'),
      el('p', { class: 'small' }, 'JLB, JSB, JSA and JLA drawn by CPJSdraw from each record’s annotated IR. It gives SVG, PNG and a table of the junctions.'),
      el('div', { class: 'row' }, records.map((record) => el('label', { class: 'check' }, el('input', { type: 'checkbox', onchange: (event) => {
        junction.inputs = event.target.checked ? [...junction.inputs, record.path] : junction.inputs.filter((path) => path !== record.path);
        refreshJunction();
      } }), record.name))),
      el('div', { class: 'field-grid' },
        textField('Where the data came from', '', (value) => { junction.label = value; refreshJunction(); }),
        pathField('CPJSdraw.pl, when it is not on PATH', '', (value) => { junction.cpjsdraw = value; refreshJunction(); })),
      more('Drawing options',
        textField('Base name of the files', '', (value) => { junction.tag = value; refreshJunction(); }, 'Default: ir_junction'),
        textField('PNG width, in pixels', '', (value) => { junction.width = value; refreshJunction(); }, 'Default: 3200'),
        pathField('A tree that sets the species order', '', (value) => { junction.tree = value; refreshJunction(); }),
        textField('The species on that tree the others are ordered from', '', (value) => { junction.order_from = value; refreshJunction(); }),
        checkField('Leave pseudogenes out of the figure', false, (value) => { junction.drop_pseudo = value; refreshJunction(); })),
      junctionBox.line, junctionBox.start));
  drawBox.start.textContent = 'Draw';
  junctionBox.start.textContent = 'Draw the junctions';
  refreshDraw();
  refreshJunction();
}

/* Annotate --------------------------------------------------------------------------------------- */

async function annotatePage(query) {
  if (!session.modules) session.modules = await get('/api/modules');
  await loadTools();
  const module = session.modules.modules.find((item) => item.id === 'plastid_annotate');
  const form = { inputs: query.get('input') ? [query.get('input')] : [], reference_set: 'references', exclude: [], label: '', dataset: '', docx: false, no_cache: false, sets: {} };
  const box = commandBox('annotate', () => form);
  const refresh = () => box.refresh();
  page('annotate',
    el('header', null, el('h1', null, 'Annotate'), el('p', { class: 'lead' }, 'Transfer genes to a sequence from GenBank records you choose, same family first. The toolkit ships no plastid records.')),
    el('div', { class: 'columns' },
      el('div', { class: 'grow' },
        el('section', { class: 'card' },
          pathField('Sequence', form.inputs[0] || '', (value) => { form.inputs = value ? [value] : []; refresh(); }, { start: 'records' }),
          textField('Where the data came from', '', (value) => { form.label = value; refresh(); }, 'Printed with every table.'),
          textField('Dataset name', '', (value) => { form.dataset = value; refresh(); }, 'The folder inside runs.')),
        el('section', { class: 'card' }, el('div', { class: 'row between' }, el('h2', null, 'Reference set'), el('a', { href: '#/references' }, 'Manage references')),
          pathField('Folder of GenBank records', form.reference_set, (value) => { form.reference_set = value; refresh(); }, { start: 'references', folders: true }),
          listField('Leave out of the search', (value) => { form.exclude = value; refresh(); }, 'Accessions, comma-separated. The run names what it left out.')),
        el('section', { class: 'card' }, el('h2', null, 'Parameters'),
          el('div', { class: 'row' }, session.tools.extras.docx ? checkField('Word report', form.docx, (value) => { form.docx = value; refresh(); }) : null,
            checkField('Reuse stored results', true, (value) => { form.no_cache = !value; refresh(); })),
          el('div', { class: 'field-grid' }, module.params.filter((param) => !param.flag && param.type !== 'dir').map((param) => paramField(param, form.sets, refresh))))),
      el('aside', { class: 'card side' }, el('h2', null, 'The plan'),
        el('p', { class: 'small' }, 'The reference set is audited first: a blocked record stops the run and is named. A record that already carries genes is refused unless you ask to replace them.'),
        box.start, el('div', { class: 'small' }, el('b', null, 'The same run on the command line')), box.line)));
  box.start.textContent = 'Annotate';
  refresh();
}

/* Assemble --------------------------------------------------------------------------------------- */

async function assemblePage() {
  if (!session.modules) session.modules = await get('/api/modules');
  const step = session.modules.assembly;
  const form = { kind: 'short', where: 'here', reads_1: '', reads_2: '', hifi: '', sample: '', label: '', prepared: '', result: '', source_name: '', retrieved_at: '', source_url: '', compare_to: '', skip_resource_check: false, sets: {} };
  const command = () => (form.where === 'here' ? 'assemble run' : (form.where === 'prepare' ? 'assemble prepare' : 'assemble import'));
  const fields = el('section', { class: 'card' });
  const notes = el('div');
  const area = el('aside', { class: 'card side' });
  const redraw = () => {
    const box = commandBox(command(), () => (form.kind === 'hifi' ? { ...form, reads_1: '', reads_2: '' } : { ...form, hifi: '' }));
    const refresh = () => box.refresh();
    const reads = form.kind === 'hifi'
      ? [pathField('HiFi reads', form.hifi, (value) => { form.hifi = value; refresh(); }, { start: 'reads' })]
      : [pathField('First reads', form.reads_1, (value) => { form.reads_1 = value; refresh(); }, { start: 'reads' }),
        pathField('Second reads', form.reads_2, (value) => { form.reads_2 = value; refresh(); }, { start: 'reads' })];
    fields.replaceChildren(el('h2', null, form.where === 'import' ? 'The result from elsewhere' : 'Reads'), ...(form.where === 'import'
      ? [pathField('prepared.json', form.prepared, (value) => { form.prepared = value; refresh(); }),
        pathField('Folder of the assembler’s output', form.result, (value) => { form.result = value; refresh(); }, { folders: true }),
        textField('Where it ran', form.source_name, (value) => { form.source_name = value; refresh(); }, 'For example usegalaxy.eu.'),
        textField('Date retrieved', form.retrieved_at, (value) => { form.retrieved_at = value; refresh(); }, 'YYYY-MM-DD'),
        textField('Link to the run (optional)', form.source_url, (value) => { form.source_url = value; refresh(); })]
      : [...reads, textField('Sample name', form.sample, (value) => { form.sample = value; refresh(); })]),
    textField('Where the reads came from', form.label, (value) => { form.label = value; refresh(); }),
    form.where === 'prepare' ? null : pathField('A sequence to compare with (optional)', form.compare_to, (value) => { form.compare_to = value; refresh(); }, { start: 'records' }),
    more('Step parameters', step.params.filter((param) => param.name !== 'data_label').map((param) => paramField(param, form.sets, refresh))),
    form.where === 'here' ? checkField('Run even if this computer looks too small', form.skip_resource_check, (value) => { form.skip_resource_check = value; refresh(); }) : null);
    const needed = form.kind === 'hifi' ? ['oatk', 'minimap2'] : ['getorganelle', 'fastp', 'minimap2'];
    const missing = session.tools.tools.filter((tool) => needed.includes(tool.name) && tool.status === 'missing').map((tool) => tool.name);
    notes.replaceChildren(...(form.where === 'here' && missing.length ? [el('div', { class: 'note-warn' },
      `${missing.join(', ')} not found. Install them, or run the same assembly elsewhere.`,
      el('button', { class: 'btn way', type: 'button', onclick: () => { form.where = 'prepare'; draw(); } }, 'Prepare to run it elsewhere'))] : []));
    const machine = session.tools.machine;
    const needs = (name) => step.params.find((param) => param.name === name).default;
    const fits = machine.ram_gb >= Number(form.sets.min_ram_gb ?? needs('min_ram_gb'));
    area.replaceChildren(el('h2', null, 'The plan'),
      form.where === 'here' ? el('ol', { class: 'plan' },
        el('li', null, el('b', { class: `word-${fits ? 'ok' : 'warn'}` }, 'Memory'), el('span', null, `${machine.ram_gb.toFixed(1)} GB on this computer`, el('small', null, `The step stops below ${form.sets.min_ram_gb ?? needs('min_ram_gb')} GB.`))),
        el('li', null, el('b', { class: 'word-quiet' }, 'Disk'), el('span', null, `${machine.disk_free_gb.toFixed(1)} GB free`, el('small', null, `The step needs ${form.sets.disk_factor ?? needs('disk_factor')} times the size of the reads, and checks that itself.`))),
        el('li', null, el('b', { class: 'word-quiet' }, 'Cores'), el('span', null, `${machine.cores} here`, el('small', null, `The tools use ${form.sets.threads ?? needs('threads')}.`)))) : null,
      el('p', { class: 'small' }, form.where === 'here' ? 'The reads are mapped back and every join of the circle is read. A sequence moves on only if the reads cross every join.'
        : (form.where === 'prepare' ? 'Quality check and hashes. Nothing leaves this computer; you run the same tool on the service you choose.'
          : 'The result is checked against the same reads, and the run says where it ran.')),
      box.start, notes, el('div', { class: 'small' }, el('b', null, 'The same step on the command line')), box.line);
    box.start.textContent = { here: 'Assemble', prepare: 'Prepare the reads', import: 'Import the result' }[form.where];
    box.refresh(form.where === 'here' && missing.length ? { words: 'tools' } : null);
  };
  await loadTools();
  const choices = el('div', { class: 'row' });
  const draw = () => {
    const seg = (key, options) => el('div', { class: 'seg', role: 'group' }, options.map(([id, words]) => el('button', { type: 'button', 'aria-pressed': String(form[key] === id), onclick: () => { form[key] = id; draw(); } }, words)));
    choices.replaceChildren(seg('kind', [['short', 'Paired short reads'], ['hifi', 'HiFi reads']]), seg('where', [['here', 'On this computer'], ['prepare', 'Prepare for elsewhere'], ['import', 'Import a result']]));
    redraw();
  };
  page('assemble',
    el('header', null, el('h1', null, 'Assemble'), el('p', { class: 'lead' }, 'From reads to a checked plastome. The toolkit sends nothing anywhere: a run elsewhere is one you start yourself.')),
    choices, el('div', { class: 'columns' }, el('div', { class: 'grow' }, fields), area));
  draw();
}

/* References ------------------------------------------------------------------------------------- */

async function referencesPage() {
  if (!session.modules) session.modules = await get('/api/modules');
  const naming = session.modules.modules.find((module) => module.id === 'gene_normalize');
  const form = { taxon: '', max_records: '', index: '', target: '', mode: 'tiers', level: 'family', per_tier: '', accessions: '', groups: '', files: [], note: '', to: 'references', reference_set: 'references', family: '', rule: '', exceptions: '', no_confirm: false, confirm_exception: '', confirmed_by: '', reason: '' };
  const pin = () => ({ ...form, accessions: form.accessions.split(',').map((item) => item.trim()).filter(Boolean) });
  const audit = () => ({ ...form, confirm_exception: '', confirmed_by: '', reason: '' });
  const editing = { records: '', reference_file: '', to: 'references/editing_reference.json' };

  const standard = { star: '', records: [], label: '', dataset: '', sets: {} };
  const starBox = commandBox('run', () => ({
    modules: [naming.id], inputs: standard.records.map((record) => record.path).filter((path) => path !== standard.star),
    label: standard.label, dataset: standard.dataset, sets: { ...standard.sets, [naming.reference_param]: standard.star },
  }));
  const refreshStar = () => starBox.refresh(!standard.star);
  const list = el('div', { class: 'picks' });
  const readSet = async () => {
    list.replaceChildren(el('p', { class: 'small' }, 'Reading the set…'));
    const { records } = await get('/api/records', { folder: form.reference_set });
    standard.records = await Promise.all(records.map((entry) => get('/api/record', { path: entry.path })));
    standard.star = '';
    list.replaceChildren(...(standard.records.length ? standard.records.map((record) => el('label', { class: 'check' },
      el('input', { type: 'radio', name: 'star', onchange: () => { standard.star = record.path; refreshStar(); } }),
      record.unreadable ? `${record.path}: ${record.unreadable}` : el('span', null, species(record.organism), ' ', el('code', { class: 'small' }, record.accession || ''), ` · ${record.genes.total} genes`)))
      : [el('p', { class: 'small' }, 'No record in that folder.')]));
    refreshStar();
  };

  page('references',
    el('header', null, el('h1', null, 'References'), el('p', { class: 'lead' }, 'Find records, read what the audit reports, and choose one as the naming standard. That choice is yours; the toolkit never makes it.')),
    el('section', { class: 'card' }, el('h2', null, '1. Find records'),
      el('div', { class: 'note-plain' }, 'This is the one place the toolkit goes online, and only when you press a button that says so. It reads public lists and records from NCBI and sends none of your data.'),
      el('h3', null, 'The list for a taxon'),
      el('div', { class: 'field-grid' },
        textField('Taxon', '', (value) => { form.taxon = value; }, 'For example Iridaceae.'),
        textField('Most records to ask for', '', (value) => { form.max_records = value; }, 'Default: 200'),
        textField('Where the list is kept', '', (value) => { form.index = value; }, 'Default: .panel/index.json')),
      el('div', { class: 'row' }, jobButton('panel refresh', () => form, 'Ask NCBI for the list')),
      el('h3', null, 'Candidates for your genome'),
      el('div', { class: 'field-grid' },
        textField('Your target genome', '', (value) => { form.target = value; }, 'An accession on that list.'),
        selectField('How to choose', ['tiers', 'all'], form.mode, (value) => { form.mode = value; }),
        selectField('How far to look, with all', ['family', 'species', 'genus', 'order'], form.level, (value) => { form.level = value; }),
        textField('How many from each tier, with tiers', '', (value) => { form.per_tier = value; }, 'Default: 3')),
      el('div', { class: 'row' }, jobButton('panel propose', () => form, 'Propose candidates')),
      el('h3', null, 'Fetch and pin'),
      el('div', { class: 'field-grid' },
        textField('Accessions to fetch', '', (value) => { form.accessions = value; }, 'Comma-separated.'),
        pathField('Or a file of sample groups', '', (value) => { form.groups = value; }),
        pathField('Where the records go', form.to, (value) => { form.to = value; }, { start: 'references', folders: true })),
      el('div', { class: 'row' }, jobButton('panel pin', pin, 'Fetch and pin these')),
      more('Pin a record made elsewhere',
        pathField('The record', '', (value) => { form.files = value ? [value] : []; }),
        textField('Which service made it, and when', '', (value) => { form.note = value; }, 'For example PMGA web, 2026-10-03.'),
        jobButton('panel pin', pin, 'Pin this record'))),
    el('section', { class: 'card' }, el('h2', null, '2. Audit a set'),
      pathField('Folder of GenBank records', form.reference_set, (value) => { form.reference_set = value; }, { start: 'references', folders: true }),
      el('p', { class: 'small' }, 'The audit blocks what cannot be used and reports what a person must know. Differences are reported, not repaired.'),
      el('div', { class: 'row' }, jobButton('reference audit', audit, 'Audit the set', 'btn go'), jobButton('reference verify', () => form, 'Check the files')),
      more('Audit options',
        pathField('The file of pinned exceptions', '', (value) => { form.exceptions = value; }),
        checkField('Skip the second reading; the report then says not measured', false, (value) => { form.no_confirm = value; })),
      more('Accept one reported finding',
        textField('File and check', '', (value) => { form.confirm_exception = value; }, 'FILE:CHECK, as the report writes it. Editing the file voids it.'),
        textField('Who confirms it', '', (value) => { form.confirmed_by = value; }),
        textField('Why it is accepted', '', (value) => { form.reason = value; }),
        jobButton('reference audit', () => form, 'Pin this exception')),
      more('Which record may annotate which',
        textField('The family the set was chosen at', '', (value) => { form.family = value; }, 'For example Iridaceae.'),
        textField('The sentence it was chosen by', '', (value) => { form.rule = value; }, 'Default: same family first'),
        pathField('Where to write it', '', (value) => { form.to_set = value; }),
        jobButton('reference leave-one-out', () => ({ ...form, to: form.to_set }), 'Write it, each record left out of its own set'))),
    el('section', { class: 'card' }, el('h2', null, '3. Choose the naming standard'),
      el('p', { class: 'small' }, 'Every other record of the set is named after the one you mark here. Unifying writes a new set with its own checksum. Read the audit first: the record you choose is audited again, and one that cannot be used stops the run.'),
      el('div', { class: 'row' }, el('button', { class: 'btn', type: 'button', onclick: () => readSet().catch((error) => list.replaceChildren(fail(error))) }, 'Read the set above')),
      list,
      el('div', { class: 'field-grid' },
        textField('Where the data came from', '', (value) => { standard.label = value; refreshStar(); }, 'Printed with every table.'),
        textField('Dataset name', '', (value) => { standard.dataset = value; refreshStar(); }, 'The folder inside runs.'),
        naming.params.filter((param) => ['mode', 'on_conflict'].includes(param.name)).map((param) => paramField(param, standard.sets, refreshStar))),
      starBox.line, starBox.start),
    el('section', { class: 'card' }, el('h2', null, 'The reference behind RNA-editing predictions'),
      el('p', { class: 'small' }, 'A rebuild of the PREP-Cp training reference from seven pinned records. It is not the original PREP-Cp reference, and the results say so.'),
      el('div', { class: 'row' },
        jobButton('editing-reference verify', () => ({ reference_file: editing.reference_file, records: editing.records }), 'Verify it', 'btn go'),
        jobButton('editing-reference build', () => editing, 'Rebuild it')),
      more('Which files',
        pathField('Folder of the seven pinned records', '', (value) => { editing.records = value; }, { folders: true }),
        pathField('A reference file to verify', '', (value) => { editing.reference_file = value; }),
        textField('Where a rebuild is written', editing.to, (value) => { editing.to = value; }))));
  starBox.start.textContent = 'Unify the names';
  refreshStar();
}

/* Tools ------------------------------------------------------------------------------------------ */

async function toolsPage() {
  const body = el('tbody');
  const redraw = (tools) => body.replaceChildren(...tools.tools.flatMap((tool) => {
    const fine = tool.status === 'ok';
    return [el('tr', null,
      el('td', { class: 'num' }, tool.name),
      el('td', null, el('span', { class: `word-${fine ? 'ok' : 'warn'}` }, { ok: 'Found', missing: 'Missing', other_version: 'Another version' }[tool.status] || tool.status)),
      el('td', { class: 'num' }, tool.version || 'none'), el('td', { class: 'num' }, tool.pinned || ''),
      el('td', null, tool.used_by.join(', '))),
    fine ? null : el('tr', null, el('td', { colspan: 5, class: 'small' }, tool.message))].filter(Boolean);
  }));
  const tools = await loadTools();
  page('tools',
    el('header', { class: 'page-head' },
      el('div', null, el('h1', null, 'Tools'), el('p', { class: 'lead' }, 'Version is part of a result, so each tool is checked against the version the toolkit’s numbers were measured with. A different build is reported, never used quietly.')),
      el('button', { class: 'btn', type: 'button', onclick: async () => { session.tools = await post('/api/tools/refresh'); toolsPage(); } }, 'Check again')),
    el('p', { class: 'small' }, `This computer: ${tools.machine.ram_gb.toFixed(1)} GB memory, ${tools.machine.cores} cores, ${tools.machine.disk_free_gb.toFixed(1)} GB free in your folder.`),
    el('p', { class: 'small' }, 'A tool kept somewhere else can be named on the analysis that uses it: its parameters have a field for the tool. “If asked” means only a chosen option needs that tool.'),
    el('div', { class: 'table-box' }, el('table', null,
      el('thead', null, el('tr', null, ['Tool', 'Status', 'Found', 'Pinned', 'Needed by'].map((name) => el('th', { scope: 'col' }, name)))), body)));
  redraw(tools);
}

/* Addresses -------------------------------------------------------------------------------------- */

async function route() {
  const [path, search] = location.hash.replace(/^#\/?/, '').split('?');
  const parts = path.split('/').filter(Boolean).map(decodeURIComponent);
  const query = new URLSearchParams(search || '');
  try {
    if (!session.about) session.about = await get('/api/app');
    const [name, ...rest] = parts;
    session.active = { job: 'runs', run: 'runs' }[name] || (session.about.pages.includes(name) ? name : 'start');
    if (name === 'records') await recordsPage();
    else if (name === 'analyse') await analysePage();
    else if (name === 'job') await jobPage(rest[0]);
    else if (name === 'runs') await runsPage();
    else if (name === 'run') await runPage(rest[0], rest[1], rest[2]);
    else if (name === 'export') await exportPage(query);
    else if (name === 'annotate') await annotatePage(query);
    else if (name === 'assemble') await assemblePage();
    else if (name === 'references') await referencesPage();
    else if (name === 'tools') await toolsPage();
    else await startPage();
  } catch (error) {
    page('', el('h1', null, 'That did not work'), fail(error));
  }
}

window.addEventListener('hashchange', route);
route();
