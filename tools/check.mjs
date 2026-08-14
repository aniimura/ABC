#!/usr/bin/env node
/**
 * check.mjs — 本プロジェクト唯一のゲート
 *
 * 検査するもの:
 *   [A] 1_Structured の構造化文書  … S1-S6(規約は 1_Structured/README.md §4)
 *   [B] Lean プロジェクト          … ビルド / sorry 台帳 / axiom 禁止(G4)
 *
 * 使い方:
 *   node tools/check.mjs                 全部
 *   node tools/check.mjs --structured    [A] のみ
 *   node tools/check.mjs --lean          [B] のみ
 *   node tools/check.mjs --selftest      器具自身の較正
 *
 * 終了コード: 0 = 全 PASS / 1 = NG あり
 *
 * ★原則: 器具は自分の欠陥を報告できなければならない。
 *   非ゲート指標(件数等)も必ず数値を印字し、「緑1語」で隠さない。
 */

import {
  readFileSync, readdirSync, statSync, existsSync,
  mkdtempSync, mkdirSync, writeFileSync, rmSync,
} from 'node:fs';
import { execFileSync } from 'node:child_process';
import { join, dirname, relative, basename } from 'node:path';
import { fileURLToPath } from 'node:url';
import { tmpdir } from 'node:os';

const ROOT = dirname(dirname(fileURLToPath(import.meta.url)));
const SOURCE_DIR = join(ROOT, 'ResearchPaper', '0_Source');
const STRUCT_DIR = join(ROOT, 'ResearchPaper', '1_Structured');
const PAPERS_JSON = join(ROOT, 'ResearchPaper', 'papers.json');
const LEAN_DIR = join(ROOT, 'lean');
const LEAN_SRC = join(LEAN_DIR, 'ABC3');

/** 較正デモだけは axiom を持つことを許す(それが実演の内容だから) */
const AXIOM_EXEMPT = ['Meta/Calibration.lean']; // lean/ABC3 からの相対

let NG = 0;
const ng = (where, msg) => { NG++; console.log(`  NG  ${where}\n      ${msg}`); };
const ok = (msg) => console.log(`  ok  ${msg}`);
const h1 = (s) => console.log(`\n=== ${s} ===`);

// ────────────────────────────────────────────────────────────────
// テキスト正規化 / HTML の最小パーサ
// ────────────────────────────────────────────────────────────────

const ENTITIES = {
  amp: '&', lt: '<', gt: '>', quot: '"', apos: "'", nbsp: ' ',
  ldquo: '“', rdquo: '”', lsquo: '‘', rsquo: '’',
  minus: '−', times: '×', ge: '≥', le: '≤',
  rarr: '→', larr: '←', harr: '↔',
  sube: '⊆', supe: '⊇', sub: '⊂', sup: '⊃',
  cong: '≅', or: '∨', and: '∧',
  Gamma: 'Γ', alpha: 'α', chi: 'χ', prime: '′',
  sect: '§', hellip: '…', mdash: '—', ndash: '–',
};

function decodeEntities(s) {
  return s
    .replace(/&#x([0-9a-fA-F]+);/g, (_, h) => String.fromCodePoint(parseInt(h, 16)))
    .replace(/&#(\d+);/g, (_, d) => String.fromCodePoint(parseInt(d, 10)))
    .replace(/&([a-zA-Z]+);/g, (m, name) => (name in ENTITIES ? ENTITIES[name] : m));
}

/**
 * 照合用の正規化。両側に同じものを適用する。
 *
 * ここで吸収するのは **レンダラ由来の差** だけ——合字・引用符・ダッシュ・
 * 記号の分解(pdftotext は ≅ を "∼=" や "=∼" として出力する)。
 * **論文由来の差**(装飾の有無)はここで吸収してはならない——それは
 * `data-txt` で明示するもので、自動で潰すと事実2の検出力が消える。
 */
function normalize(s) {
  return decodeEntities(s)
    .replace(/ﬀ/g, 'ff').replace(/ﬁ/g, 'fi').replace(/ﬂ/g, 'fl')
    .replace(/ﬃ/g, 'ffi').replace(/ﬄ/g, 'ffl')
    .replace(/­/g, '')
    .replace(/[‐-―]/g, '-')
    .replace(/[‘’]/g, "'").replace(/[“”]/g, '"')
    .replace(/∼=|=∼/g, '≅')
    .replace(/\s+/g, ' ')
    .trim();
}

/**
 * 照合は空白を完全に無視して行う。
 *
 * 理由: pdftotext は上付き・下付きを別のテキスト run として出力するため、
 * 添字の直後に空白が入る(原文 "Γ_K." → 出力 "ΓK ." / "Γ_K-module" →
 * "ΓK -module")。この空白は組版の産物であって原文の一部ではないので、
 * 照合の対象にしない。空白を無視しても、数百文字の部分文字列一致という
 * 検査の識別力はほとんど落ちない。
 */
const squash = (s) => normalize(s).replace(/\s+/g, '');

/**
 * .verbatim の内側から「照合射影」を作る。
 * data-txt を持つ要素はその値に置換する(空文字なら削除)。それ以外はタグを剥がす。
 */
function matchProjection(html) {
  let s = html;
  const withTxt = /<(\w+)\b[^>]*\bdata-txt="([^"]*)"[^>]*>([\s\S]*?)<\/\1>/;
  for (let guard = 0; guard < 500 && withTxt.test(s); guard++) {
    s = s.replace(withTxt, (_m, _tag, txt) => txt);
  }
  s = s.replace(/<[^>]+>/g, '');
  return squash(s);
}

/** 極小の section 抽出。属性は key="value" 形式のみを想定する。 */
function extractSections(html) {
  const out = [];
  const re = /<section\b([^>]*)>([\s\S]*?)<\/section>/g;
  let m;
  while ((m = re.exec(html)) !== null) {
    const attrs = {};
    const ar = /([\w-]+)="([^"]*)"/g;
    let a;
    while ((a = ar.exec(m[1])) !== null) attrs[a[1]] = a[2];
    const vm = /<div class="verbatim">([\s\S]*?)<\/div>/.exec(m[2]);
    out.push({ attrs, verbatimHtml: vm ? vm[1] : null });
  }
  return out;
}

function walk(dir, acc = []) {
  if (!existsSync(dir)) return acc;
  for (const name of readdirSync(dir)) {
    const p = join(dir, name);
    if (statSync(p).isDirectory()) walk(p, acc);
    else acc.push(p);
  }
  return acc;
}

// ────────────────────────────────────────────────────────────────
// PDF ページのテキスト(キャッシュ付き)
// ────────────────────────────────────────────────────────────────

/**
 * pdftotext の抽出モード。
 *
 * 既定モードは、2段組でない本文でも **行の順序を入れ替える**ことがある
 * (実測: pGC 物理 p.3 の Corollary 1.3 が "Corollary 1.3: from ΓK ." の後に
 * 本体行が来る)。`-layout` は組版上の行順を保つのでこれを回避できる。
 * 逆に `-layout` が不利な版面もありうるため、3モードを順に試し、
 * **どのモードで一致したかを必ず印字する**(既定以外での一致は、その版面が
 * 素直でないことの情報になる)。
 */
const PDF_MODES = [
  ['layout', ['-layout']],
  ['default', []],
  ['raw', ['-raw']],
];

const pageCache = new Map();
/** @returns {Array<[string, string]>|null} [モード名, squash済みテキスト] の配列 */
function pdfPageTexts(pdfPath, page) {
  const key = `${pdfPath}#${page}`;
  if (pageCache.has(key)) return pageCache.get(key);
  const out = [];
  for (const [name, flags] of PDF_MODES) {
    try {
      out.push([name, squash(execFileSync('pdftotext',
        ['-enc', 'UTF-8', ...flags, '-f', String(page), '-l', String(page), pdfPath, '-'],
        { encoding: 'utf8', maxBuffer: 32 * 1024 * 1024 }))]);
    } catch { /* このモードは使えない */ }
  }
  const val = out.length ? out : null;
  pageCache.set(key, val);
  return val;
}

// ────────────────────────────────────────────────────────────────
// [A] 1_Structured の検査
// ────────────────────────────────────────────────────────────────

const DATE_RE = /^\d{4}-\d{2}-\d{2}$/;

function checkStructured({ files = null, papersPath = PAPERS_JSON, quiet = false } = {}) {
  const before = NG;
  const reg = JSON.parse(readFileSync(papersPath, 'utf8')).papers;

  const targets = files ?? walk(STRUCT_DIR)
    .filter((p) => p.endsWith('.html'))
    .filter((p) => !basename(p).includes('.legacy.'));

  let nSection = 0;
  let nVerbatim = 0;
  const nonPrimaryMode = [];
  const nLegacy = files ? 0 : walk(STRUCT_DIR).filter((p) => basename(p).includes('.legacy.')).length;

  for (const file of targets) {
    const rel = relative(ROOT, file);
    const html = readFileSync(file, 'utf8');
    const seenIds = new Set();

    for (const { attrs, verbatimHtml } of extractSections(html)) {
      if (!/\bstatement\b/.test(attrs['class'] ?? '')) continue;
      nSection++;
      const id = attrs['id'] ?? '(id無し)';
      const at = `${rel} #${id}`;

      const required = ['id', 'data-paper', 'data-pdf-page', 'data-item', 'data-notation-checked'];
      const missing = required.filter((k) => !(k in attrs));
      if (missing.length) { ng(at, `S1 必須属性が欠落: ${missing.join(', ')}`); continue; }

      if (seenIds.has(id)) ng(at, 'S6 id が重複している');
      seenIds.add(id);

      const tag = attrs['data-paper'];
      const paper = reg[tag];
      if (!paper) { ng(at, `S2 data-paper="${tag}" が papers.json に無い`); continue; }

      const page = Number(attrs['data-pdf-page']);
      if (!Number.isInteger(page) || page < 1 || page > paper.pdfPages) {
        ng(at, `S3 data-pdf-page=${attrs['data-pdf-page']} が範囲外(1..${paper.pdfPages})`);
        continue;
      }

      const nc = attrs['data-notation-checked'];
      if (!(DATE_RE.test(nc) || nc === 'none')) {
        ng(at, `S5 data-notation-checked="${nc}" は YYYY-MM-DD か "none" でなければならない`);
      }

      if (verbatimHtml === null) { ng(at, 'S4 .verbatim が無い'); continue; }
      nVerbatim++;

      const pdfPath = join(SOURCE_DIR, `${paper.file}.pdf`);
      if (!existsSync(pdfPath)) { ng(at, `S4 PDF が見つからない: ${pdfPath}`); continue; }

      const texts = pdfPageTexts(pdfPath, page);
      if (texts === null) { ng(at, `S4 pdftotext が失敗した(物理 p.${page})`); continue; }

      const proj = matchProjection(verbatimHtml);
      const hit = texts.find(([, t]) => t.includes(proj));
      if (hit) {
        if (hit[0] !== PDF_MODES[0][0]) nonPrimaryMode.push(`${at} (${hit[0]})`);
      } else {
        // 最も長く一致したモードで、どこまで一致したかを二分探索で示す
        let best = { mode: '', lo: -1 };
        for (const [mode, t] of texts) {
          let lo = 0;
          let hi = proj.length;
          while (lo < hi) {
            const mid = Math.ceil((lo + hi) / 2);
            if (t.includes(proj.slice(0, mid))) lo = mid; else hi = mid - 1;
          }
          if (lo > best.lo) best = { mode, lo };
        }
        ng(at,
          `S4 逐語が物理 p.${page} に見つからない(最良 ${best.mode} モードで先頭 ${best.lo}/${proj.length} 文字まで一致)\n` +
          `      一致した末尾: ...${JSON.stringify(proj.slice(Math.max(0, best.lo - 60), best.lo))}\n` +
          `      次に来るはず: ${JSON.stringify(proj.slice(best.lo, best.lo + 60))}`);
      }
    }
  }

  if (!quiet) {
    if (nonPrimaryMode.length) {
      console.log(`  -- layout 以外のモードで一致した逐語: ${nonPrimaryMode.length} 件`);
      for (const x of nonPrimaryMode) console.log(`     ${x}`);
    }
    console.log(`  -- 構造単位 ${nSection} 件 / 逐語照合 ${nVerbatim} 件 / 対象文書 ${targets.length} 件` +
      (files ? '' : ` / legacy 除外 ${nLegacy} 件`));
    if (NG === before) ok('1_Structured: S1-S6 すべて PASS');
  }
  return NG - before;
}

// ────────────────────────────────────────────────────────────────
// [B-2] Lean ソースの台帳検査(G1-Lean / G2 / G4 / Gap)
//
// ★分業: **型が付くかどうかは Lean が検査する**。ここが見るのは
//   「必要な宣言が存在するか」という帳簿だけ。
//   `Foo.nonvacuous` が本当に `Nonempty Foo` を証明しているかは `lake build`
//   が保証する——check.mjs はその存在しか見ていない。この限界を明示する。
// ────────────────────────────────────────────────────────────────

const DECL_RE = /^\s*(?:@\[[^\]]*\]\s*)?(?:private\s+|protected\s+|noncomputable\s+)*(theorem|lemma|def|abbrev|instance|structure|inductive|example)\s+([A-Za-z_][\w.'!?₀-₉]*)/;

/** Lean ソース木を走査して、宣言名の集合と structure 一覧を作る。 */
function scanLeanTree(dir) {
  const files = walk(dir).filter((p) => p.endsWith('.lean'));
  const names = new Set();
  const structures = [];   // { name, file, line, bucket }
  const decls = [];        // { kind, name, file, line, bucket }
  const axioms = [];       // { file, line }
  const texts = new Map();

  for (const f of files) {
    const src = readFileSync(f, 'utf8');
    texts.set(f, src);
    const relf = relative(dir, f).replace(/\\/g, '/');
    const bucket = relf.split('/')[0];               // Interface / Skeleton / Found / Gap / Meta
    // コメントを潰してから宣言を拾う(docstring 内の例を宣言と誤認しない)
    const stripped = src
      .replace(/\/-[\s\S]*?-\//g, (m) => m.replace(/[^\n]/g, ' '))
      .replace(/--[^\n]*/g, (m) => ' '.repeat(m.length));
    stripped.split('\n').forEach((line, i) => {
      const m = DECL_RE.exec(line);
      if (m) {
        names.add(m[2]);
        decls.push({ kind: m[1], name: m[2], file: f, line: i + 1, bucket });
        if (m[1] === 'structure' || m[1] === 'inductive') {
          structures.push({ name: m[2], file: f, line: i + 1, bucket });
        }
      }
      if (/^\s*axiom\s+/.test(line)) axioms.push({ file: f, line: i + 1 });
    });
  }
  return { files, names, structures, decls, axioms, texts };
}

/** 1_Structured から「論文タグ → section id の集合」を集める。 */
function collectStructuredIds() {
  const map = new Map();
  const files = walk(STRUCT_DIR)
    .filter((p) => p.endsWith('.html'))
    .filter((p) => !basename(p).includes('.legacy.'));
  for (const f of files) {
    for (const { attrs } of extractSections(readFileSync(f, 'utf8'))) {
      if (!/\bstatement\b/.test(attrs['class'] ?? '')) continue;
      const paper = attrs['data-paper'];
      const id = attrs['id'];
      if (!paper || !id) continue;
      if (!map.has(paper)) map.set(paper, new Set());
      map.get(paper).add(id);
    }
  }
  return map;
}

/**
 * @param {object} o
 * @param {string} o.dir            走査するルート(本番は lean/ABC3、selftest は fixture)
 * @param {string[]} o.axiomExempt  bucket 相対パスで axiom を許すファイル
 */
function checkLeanLedger({ dir, axiomExempt = [], papersPath = PAPERS_JSON, quiet = false } = {}) {
  const before = NG;
  if (!existsSync(dir)) { ng(relative(ROOT, dir), 'ディレクトリが無い'); return NG - before; }

  const reg = JSON.parse(readFileSync(papersPath, 'utf8')).papers;
  const structuredIds = collectStructuredIds();
  const { names, structures, decls, axioms, texts } = scanLeanTree(dir);
  const at = (x) => `${relative(ROOT, x.file)}:${x.line}`;

  // ── G4: axiom 禁止
  for (const a of axioms) {
    const relf = relative(dir, a.file).replace(/\\/g, '/');
    if (!axiomExempt.includes(relf)) {
      ng(`${relative(ROOT, a.file)}:${a.line}`,
        'G4 axiom 宣言は禁止(未構築の基礎は Interface/ の structure で受ける)');
    }
  }

  // ── G2: Interface の structure は nonvacuous witness か waiting を持つ
  const waiting = [];
  for (const s of structures.filter((x) => x.bucket === 'Interface')) {
    const hasWitness = names.has(`${s.name}.nonvacuous`);
    const hasWaiting = names.has(`${s.name}.waiting`);
    if (!hasWitness && !hasWaiting) {
      ng(at(s), `G2 非空虚 witness が無い: \`${s.name}.nonvacuous : Nonempty ${s.name}\` を構成するか、` +
                `構成できないなら \`${s.name}.waiting : ABC3.Meta.WaitingFor\` で何を待っているか書く`);
    } else if (hasWaiting && !hasWitness) {
      const src = texts.get(s.file) ?? '';
      const m = new RegExp(`${s.name.replace(/[.*+?^${}()|[\\]\\\\]/g, '\\\\$&')}\\.waiting[\\s\\S]{0,400}?what\\s*:=\\s*"([^"]*)"[\\s\\S]{0,200}?trackB\\s*:=\\s*"([^"]*)"`).exec(src);
      waiting.push(m ? `${s.name} — ${m[1]} → ${m[2]}` : `${s.name} — (waiting の内容を読めなかった)`);
    }
  }

  // ── Gap: 各 structure は GapRecord を伴う(falsifier は型が必須にしている)
  const gaps = [];
  for (const s of structures.filter((x) => x.bucket === 'Gap')) {
    if (!names.has(`${s.name}.record`)) {
      ng(at(s), `Gap の記録が無い: \`${s.name}.record : ABC3.Meta.GapRecord\` を書く` +
                '(falsifier は GapRecord の必須フィールド——書けないうちは ③ ではない)');
      continue;
    }
    const src = texts.get(s.file) ?? '';
    const idx = src.indexOf(`${s.name}.record`);
    const m = /classification\s*:=\s*\.?(\w+)[\s\S]{0,400}?falsifier\s*:=\s*"([^"]*)"/
      .exec(idx >= 0 ? src.slice(idx) : src);
    gaps.push(m ? `${s.name} — ${m[1]}: ${m[2]}` : `${s.name} — (record の内容を読めなかった)`);
  }

  // ── G3: load-bearing と印を付けた宣言には負の対照が要る
  for (const d of decls.filter((x) => x.bucket === 'Skeleton')) {
    if (!names.has(`${d.name}.loadBearing`)) continue;
    if (!names.has(`${d.name}.negControl`)) {
      ng(at(d), `G3 負の対照が無い: \`${d.name}.negControl : ABC3.Meta.NegControl\` を書く` +
                '(性質を1つだけ落とした対照が破れることを確認する。破れないならその性質は効いていない)');
    }
  }

  // ── 主語の分離: Skeleton(原典の主張)は Check(我々のモデルの検査)に依存できない
  for (const [f, src] of texts) {
    const relf = relative(dir, f).replace(/\\/g, '/');
    if (relf.split('/')[0] !== 'Skeleton') continue;
    src.split('\n').forEach((line, i) => {
      if (/^\s*import\s+ABC3\.Check\b/.test(line)) {
        ng(`${relative(ROOT, f)}:${i + 1}`,
          '主語の分離: Skeleton(原典の主張)から Check(我々のモデルの検査)を import してはならない');
      }
    });
  }

  // ── G1-Lean: Skeleton の宣言は Source を伴い、その locator が実在する
  const SRC_RE = /\.src\b[\s\S]{0,400}?paper\s*:=\s*"([^"]*)"[\s\S]{0,300}?pdfPage\s*:=\s*(\d+)[\s\S]{0,300}?sectionId\s*:=\s*"([^"]*)"/;
  let nSrcOk = 0;
  for (const d of decls.filter((x) => x.bucket === 'Skeleton')) {
    // 台帳の付随宣言そのものには出典を要求しない
    if (/\.(src|nonvacuous|waiting|record|loadBearing|negControl)$/.test(d.name)) continue;
    if (!['theorem', 'lemma', 'def', 'abbrev', 'structure'].includes(d.kind)) continue;
    if (!names.has(`${d.name}.src`)) {
      ng(at(d), `G1 出典が無い: \`${d.name}.src : ABC3.Meta.Source\` を書く`);
      continue;
    }
    const src = texts.get(d.file) ?? '';
    const idx = src.indexOf(`${d.name}.src`);
    const m = SRC_RE.exec(idx >= 0 ? src.slice(idx) : src);
    if (!m) { ng(at(d), `G1 \`${d.name}.src\` の中身を読めなかった(1行の正準形で書くこと)`); continue; }
    const [, paper, pageStr, sectionId] = m;
    const page = Number(pageStr);
    const p = reg[paper];
    if (!p) { ng(at(d), `G1 paper="${paper}" が papers.json に無い`); continue; }
    if (page < 1 || page > p.pdfPages) { ng(at(d), `G1 pdfPage=${page} が範囲外(1..${p.pdfPages})`); continue; }
    const ids = structuredIds.get(paper);
    if (!ids || !ids.has(sectionId)) {
      ng(at(d), `G1 sectionId="${sectionId}" が 1_Structured(${paper})に無い` +
                `——先に構造化するか、id を直す`);
      continue;
    }
    nSrcOk++;
  }

  if (!quiet) {
    console.log(`  -- Lean 宣言 ${decls.length} / structure ${structures.length} / axiom ${axioms.length} 件` +
                `(免除 ${axiomExempt.length} ファイル)`);
    const nCheck = decls.filter((x) => x.bucket === 'Check').length;
    console.log(`  -- Skeleton の出典照合: ${nSrcOk} 件 OK / Check(我々のモデルの検査): ${nCheck} 宣言`);
    console.log(`  -- Interface 実装待ち: ${waiting.length} 件(= Track B の作業キュー)`);
    for (const w of waiting) console.log(`     ${w}`);
    console.log(`  -- Gap(飛躍): ${gaps.length} 件`);
    for (const g of gaps) console.log(`     ${g}`);
    console.log('  -- 注意: ここは宣言の**存在**しか見ていない。型の正しさは lake build が保証する');
  }
  return NG - before;
}

// ────────────────────────────────────────────────────────────────
// [B] Lean の検査
// ────────────────────────────────────────────────────────────────

function checkLean() {
  const before = NG;
  if (!existsSync(join(LEAN_DIR, 'lakefile.toml'))) { ng('lean/', 'lakefile.toml が無い'); return 1; }

  try {
    execFileSync('lake', ['build'], { cwd: LEAN_DIR, encoding: 'utf8', maxBuffer: 64 * 1024 * 1024 });
    ok('lake build 成功');
  } catch (e) {
    const out = `${e.stdout ?? ''}${e.stderr ?? ''}`;
    ng('lean/', `lake build 失敗\n${out.split('\n').slice(-15).join('\n')}`);
  }

  const leanFiles = walk(LEAN_SRC).filter((p) => p.endsWith('.lean'));

  // sorry 台帳(ゲートではない。件数を必ず印字する)
  // コメント(/- ... -/ と -- 行末)は改行を保ったまま空白に潰してから走査する。
  // ——docstring 内で "sorry" に言及しただけのものを数えると台帳が信用できなくなる。
  const stripComments = (src) => src
    .replace(/\/-[\s\S]*?-\//g, (m) => m.replace(/[^\n]/g, ' '))
    .replace(/--[^\n]*/g, (m) => ' '.repeat(m.length));

  const sorries = [];
  for (const f of leanFiles) {
    stripComments(readFileSync(f, 'utf8')).split('\n').forEach((line, i) => {
      if (/(^|[^\w.])sorry([^\w]|$)/.test(line)) sorries.push(`${relative(ROOT, f)}:${i + 1}`);
    });
  }
  console.log(`  -- sorry: ${sorries.length} 件${sorries.length ? `\n     ${sorries.join('\n     ')}` : ''}`);

  // G1-Lean / G2 / G4 / Gap の台帳検査
  checkLeanLedger({ dir: LEAN_SRC, axiomExempt: AXIOM_EXEMPT });

  if (NG === before) ok('Lean: ゲート PASS');
  return NG - before;
}

// ────────────────────────────────────────────────────────────────
// selftest — 器具自身の較正
// ────────────────────────────────────────────────────────────────

function selftest() {
  h1('selftest: 器具は壊れた入力を落とせるか');
  const tmp = mkdtempSync(join(tmpdir(), 'abc3-selftest-'));
  const papers = join(tmp, 'papers.json');
  writeFileSync(papers, readFileSync(PAPERS_JSON, 'utf8'));

  const mk = (name, body) => {
    const p = join(tmp, name);
    writeFileSync(p, `<!DOCTYPE html><html><body>${body}</body></html>`, 'utf8');
    return p;
  };
  const sec = (attrs, verbatim) =>
    `<section class="statement proposition" ${attrs}><div class="verbatim">${verbatim}</div></section>`;

  const cases = [
    ['D1 存在しないPDFページを指す locator',
      mk('d1.html', sec('id="d1" data-paper="pGC" data-pdf-page="999" data-item="X" data-notation-checked="2026-08-14"', 'whatever'))],
    ['D2 逐語が該当ページに無い',
      mk('d2.html', sec('id="d2" data-paper="pGC" data-pdf-page="3" data-item="X" data-notation-checked="2026-08-14"', 'This sentence does not appear anywhere in the paper.'))],
    ['D3 必須属性の欠落(data-notation-checked 無し)',
      mk('d3.html', sec('id="d3" data-paper="pGC" data-pdf-page="3" data-item="X"', 'x'))],
    ['D4 未登記の論文タグ',
      mk('d4.html', sec('id="d4" data-paper="NOPE" data-pdf-page="1" data-item="X" data-notation-checked="none"', 'x'))],
    ['D5 目視確認日が日付でない',
      mk('d5.html', sec('id="d5" data-paper="pGC" data-pdf-page="3" data-item="Proposition 1.1" data-notation-checked="yes"', 'The cyclotomic character'))],
  ];

  let passed = 0;
  for (const [label, file] of cases) {
    const before = NG;
    checkStructured({ files: [file], papersPath: papers, quiet: true });
    const caught = NG > before;
    NG = before; // selftest の結果は本体の NG に混ぜない
    console.log(`  ${caught ? 'ok ' : 'NG '} ${label} → ${caught ? '落とせた' : '★素通りした'}`);
    if (caught) passed++;
  }

  // D6: 正しい入力が通ること(偽陽性の検査)
  const real = join(STRUCT_DIR,
    'A Version of the Grothendieck Conjecture for p-adic Local Fields', 'section-1.html');
  let d6 = false;
  if (existsSync(real)) {
    const before = NG;
    checkStructured({ files: [real], quiet: true });
    d6 = NG === before;
    NG = before;
  }
  console.log(`  ${d6 ? 'ok ' : 'NG '} D6 正しい入力(pGC §1)は通る → ${d6 ? '通った' : '★落ちた(偽陽性)'}`);
  if (d6) passed++;

  // ── Lean 側の較正。fixture を Interface/ Gap/ に配置した木を組み立てて叩く。
  const leanCases = [
    ['D7 Interface に非空虚 witness も waiting も無い', 'Interface', 'd7-vacuous-interface.lean', true],
    ['D8 プロジェクト内の axiom 宣言', 'Interface', 'd8-axiom.lean', true],
    ['D9 Gap に GapRecord(falsifier)が無い', 'Gap', 'd9-gap-no-record.lean', true],
    ['D10 非空虚 witness を持つ Interface は通る', 'Interface', 'd10-good-witness.lean', false],
    ['D11 waiting を書いた Interface は通る', 'Interface', 'd11-waiting.lean', false],
    ['D12 load-bearing なのに負の対照が無い', 'Skeleton', 'd12-loadbearing-no-negcontrol.lean', true],
    ['D13 load-bearing + 負の対照ありは通る', 'Skeleton', 'd13-loadbearing-ok.lean', false],
    ['D14 Skeleton が Check を import している', 'Skeleton', 'd14-skeleton-imports-check.lean', true],
  ];
  const FIXTURES = join(ROOT, 'tools', 'selftest-fixtures');
  for (const [label, bucket, fixture, shouldFail] of leanCases) {
    const root = join(tmp, `lean-${fixture}`);
    mkdirSync(join(root, bucket), { recursive: true });
    writeFileSync(join(root, bucket, fixture), readFileSync(join(FIXTURES, fixture), 'utf8'), 'utf8');
    const before = NG;
    checkLeanLedger({ dir: root, quiet: true });
    const failed = NG > before;
    NG = before;
    const good = failed === shouldFail;
    console.log(`  ${good ? 'ok ' : 'NG '} ${label} → ` +
      (shouldFail
        ? (failed ? '落とせた' : '★素通りした')
        : (failed ? '★落ちた(偽陽性)' : '通った')));
    if (good) passed++;
  }

  rmSync(tmp, { recursive: true, force: true });
  const total = cases.length + 1 + leanCases.length;
  console.log(`\n  selftest: ${passed}/${total} PASS`);
  if (passed !== total) NG++;
  return passed === total;
}

// ────────────────────────────────────────────────────────────────

const args = process.argv.slice(2);
const only = (f) => args.includes(f);
const all = args.length === 0;

console.log(`Math_ABC3 check — ${new Date().toISOString().slice(0, 19).replace('T', ' ')}`);

if (all || only('--selftest')) selftest();
if (all || only('--structured')) { h1('1_Structured'); checkStructured(); }
if (all || only('--lean')) { h1('Lean'); checkLean(); }

console.log(`\n${NG === 0 ? 'PASS' : `NG ${NG} 件`}`);
process.exit(NG === 0 ? 0 : 1);
