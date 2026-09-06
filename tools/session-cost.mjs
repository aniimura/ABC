#!/usr/bin/env node
/**
 * セッションの **token 費用**を測る —— 2026-09-06、メタ第 10 回(backlog M26/M27)。
 *
 * ★動機: 進め方の改善はこれまで **速度**しか測っていなかった。token は 1 度も測っていない。
 * 本体セッションから「メタは速度を最適化しているのか token を最適化しているのか」と
 * 直接に問われたので、**測る道具を残す**。次の起動が同じスクリプトを書き直さないために。
 *
 * ★★測り方の要点(第 10 回に踏んだ穴。ここを外すと 1.6 倍ずれる):
 *   1. **`assistant` レコードは 1 回の API 呼び出しにつき複数ある**
 *      (thinking / text / tool_use が別レコード)。`usage` は同じものが繰り返し載る。
 *      実測: レコード 38,673 に対し一意な `requestId` は **23,524**。
 *      ⇒ **`requestId` で重複排除しないと token を 1.64 倍に数える。**
 *   2. **`subagent_tokens` は請求量ではない。** `<task-notification>` が報告する値は
 *      本体の `cache_creation + output` と同じ桁(1 呼び出しあたり 2.7K 対 3.8K)であり、
 *      `cache_read`(本体では全体の 99.3%)を**含まない**。
 *      ⇒ agent を `subagent_tokens` で並べると「**新しく作った量**」で並ぶだけで、
 *        費用の順にはならない。費用は**往復回数**に比例する。
 *
 * 使い方:
 *   node tools/session-cost.mjs              # 要約(既定)
 *   node tools/session-cost.mjs --tools      # 道具別の文脈占有
 *   node tools/session-cost.mjs --agents     # agent 1 体ごと(brief の長さ・CJK 率つき)
 *   node tools/session-cost.mjs --waste      # 重複呼び出しと、束ねられる連続走
 *   node tools/session-cost.mjs --all        # 全部
 *   node tools/session-cost.mjs --json
 *   node tools/session-cost.mjs --dir <path> # 記録の置き場を明示
 *   node tools/session-cost.mjs --file <a.jsonl> [--file <b.jsonl>]
 *   node tools/session-cost.mjs --selftest
 */

import { createReadStream, readdirSync, statSync, existsSync, writeFileSync, mkdtempSync, rmSync } from 'node:fs';
import { createInterface } from 'node:readline';
import { join, dirname, basename } from 'node:path';
import { fileURLToPath } from 'node:url';
import { homedir, tmpdir } from 'node:os';
import { createHash } from 'node:crypto';

const ROOT = dirname(dirname(fileURLToPath(import.meta.url)));
const argv = process.argv.slice(2);
const has = (f) => argv.includes(f);
const val = (f) => { const i = argv.indexOf(f); return i >= 0 ? argv[i + 1] : null; };
const vals = (f) => argv.flatMap((a, i) => (a === f ? [argv[i + 1]] : []));

/** 記録の置き場を当てる。worktree から呼ばれても本体の記録に当たるよう、祖先も試す。 */
function findProjectDir(explicit) {
  if (explicit) return explicit;
  const base = join(homedir(), '.claude', 'projects');
  if (!existsSync(base)) return null;
  // ★置き場の名前は絶対パスの区切りと `_` を `-` に潰したもの。
  //   `D:\Math_ABC3` → `D--Math-ABC3`(`:` `\` `_` がすべて `-`)。
  const slug = (p) => p.replace(/[\\/:_]/g, '-');
  const tries = [];
  let p = ROOT;
  for (let i = 0; i < 8; i++) { tries.push(slug(p)); const up = dirname(p); if (up === p) break; p = up; }
  for (const t of tries) { const d = join(base, t); if (existsSync(d)) return d; }
  // 当たらなければ「いちばん新しい」置き場
  const dirs = readdirSync(base).map((n) => join(base, n)).filter((d) => statSync(d).isDirectory());
  if (!dirs.length) return null;
  return dirs.sort((a, b) => statSync(b).mtimeMs - statSync(a).mtimeMs)[0];
}

/** 読み取り専用とみなす `Bash` の先頭語(束ねられる候補を数えるためだけに使う)。 */
const RO_BASH = /^\s*(cat|head|tail|sed\s+-n|ls|wc|find|grep|rg|git\s+(log|show|diff|status)|node\s+tools\/)\b/;
const RO_TOOL = new Set(['Read', 'Grep', 'Glob']);

const NOTIF_RE = /<task-id>([^<]+)<\/task-id>[\s\S]*?<subagent_tokens>(\d+)<\/subagent_tokens><tool_uses>(\d+)<\/tool_uses><duration_ms>(\d+)<\/duration_ms>/;
const CJK_RE = /[぀-ヿ㐀-鿿＀-￯]/g;

async function scan(files) {
  const st = {
    files, calls: 0, input: 0, cacheRead: 0, cacheCreate: 0, output: 0,
    ctxSum: 0, ctxList: [], byDay: new Map(),
    toolCalls: 0, batchDist: new Map(), noToolCalls: 0,
    inChars: new Map(), inN: new Map(), outChars: new Map(), outN: new Map(),
    dup: 0, dupBy: new Map(),
    roRuns: 0, roSaved: 0, roCalls: 0,
    launches: new Map(), notifs: new Map(),
  };
  const seenReq = new Set();
  const seenUse = new Set();
  const seenCall = new Set();
  const seenNotif = new Set();
  const idName = new Map();
  const toolsPerReq = new Map();
  let roRun = 0;

  for (const file of files) {
    const rl = createInterface({ input: createReadStream(file, { encoding: 'utf8' }), crlfDelay: Infinity });
    for await (const line of rl) {
      if (!line || line.length < 30) continue;
      const isA = line.includes('"type":"assistant"');
      const isU = line.includes('"type":"user"');
      const isN = line.includes('subagent_tokens');
      if (!isA && !isU && !isN) continue;
      let o; try { o = JSON.parse(line); } catch { continue; }

      // ── agent の完了通知
      if (isN && typeof o.content === 'string') {
        const m = NOTIF_RE.exec(o.content);
        if (m) {
          const k = `${m[1]}|${m[2]}|${m[4]}`;
          if (!seenNotif.has(k)) {
            seenNotif.add(k);
            const a = st.notifs.get(m[1]) || { tokens: 0, tools: 0, ms: 0, n: 0, ts: o.timestamp || '' };
            a.tokens = Math.max(a.tokens, +m[2]); a.tools = Math.max(a.tools, +m[3]);
            a.ms = Math.max(a.ms, +m[4]); a.n++;
            st.notifs.set(m[1], a);
          }
        }
      }
      // ── agent の起動(brief 本文が載る)
      const r = o.toolUseResult;
      if (r && r.agentId && typeof r.prompt === 'string' && !st.launches.has(r.agentId)) {
        st.launches.set(r.agentId, { desc: r.description || '', prompt: r.prompt, ts: o.timestamp || '' });
      }

      const c = o.message?.content;
      // ── usage(requestId で重複排除)
      if (o.type === 'assistant') {
        const id = o.requestId || o.message?.id;
        if (Array.isArray(c)) {
          const n = c.filter((b) => b.type === 'tool_use').length;
          if (n) toolsPerReq.set(id, (toolsPerReq.get(id) || 0) + n);
        }
        if (id && !seenReq.has(id)) {
          seenReq.add(id);
          const u = o.message?.usage;
          if (u) {
            st.calls++;
            st.input += u.input_tokens || 0;
            st.cacheRead += u.cache_read_input_tokens || 0;
            st.cacheCreate += u.cache_creation_input_tokens || 0;
            st.output += u.output_tokens || 0;
            const cx = (u.cache_read_input_tokens || 0) + (u.cache_creation_input_tokens || 0) + (u.input_tokens || 0);
            st.ctxSum += cx; st.ctxList.push(cx);
            const d = (o.timestamp || '').slice(0, 10);
            const x = st.byDay.get(d) || { calls: 0, read: 0, create: 0, out: 0, ctx: 0 };
            x.calls++; x.read += u.cache_read_input_tokens || 0; x.create += u.cache_creation_input_tokens || 0;
            x.out += u.output_tokens || 0; x.ctx += cx;
            st.byDay.set(d, x);
          }
        }
      }
      if (!Array.isArray(c)) continue;
      for (const b of c) {
        if (b.type === 'tool_use') {
          if (seenUse.has(b.id)) continue;
          seenUse.add(b.id);
          idName.set(b.id, b.name);
          st.toolCalls++;
          const L = JSON.stringify(b.input || {}).length;
          st.inChars.set(b.name, (st.inChars.get(b.name) || 0) + L);
          st.inN.set(b.name, (st.inN.get(b.name) || 0) + 1);
          const key = createHash('sha1').update(`${b.name} ${JSON.stringify(b.input)}`).digest('hex');
          if (seenCall.has(key)) { st.dup++; st.dupBy.set(b.name, (st.dupBy.get(b.name) || 0) + 1); }
          else seenCall.add(key);
          const ro = RO_TOOL.has(b.name) || (b.name === 'Bash' && RO_BASH.test(String(b.input?.command ?? '')));
          if (ro) { roRun++; st.roCalls++; }
          else { if (roRun >= 2) { st.roRuns++; st.roSaved += roRun - 1; } roRun = 0; }
        } else if (b.type === 'tool_result') {
          const n = idName.get(b.tool_use_id) || '?';
          const L = typeof b.content === 'string' ? b.content.length : JSON.stringify(b.content ?? '').length;
          st.outChars.set(n, (st.outChars.get(n) || 0) + L);
          st.outN.set(n, (st.outN.get(n) || 0) + 1);
        }
      }
    }
  }
  if (roRun >= 2) { st.roRuns++; st.roSaved += roRun - 1; }
  for (const n of toolsPerReq.values()) st.batchDist.set(n, (st.batchDist.get(n) || 0) + 1);
  st.noToolCalls = st.calls - toolsPerReq.size;
  st.reqWithTools = toolsPerReq.size;
  return st;
}

const M = (v) => `${(v / 1e6).toFixed(2)}M`;
const pct = (a, b) => (b ? `${((100 * a) / b).toFixed(2)}%` : '—');

function summary(st) {
  const total = st.input + st.cacheRead + st.cacheCreate + st.output;
  const ctx = st.calls ? Math.round(st.ctxSum / st.calls) : 0;
  const out = [];
  out.push(`記録 ${st.files.length} 本 / API 呼び出し(requestId 一意) ${st.calls}`);
  out.push('');
  out.push(`  cache_read      ${M(st.cacheRead).padStart(9)} token  ${pct(st.cacheRead, total)}`);
  out.push(`  cache_creation  ${M(st.cacheCreate).padStart(9)} token  ${pct(st.cacheCreate, total)}`);
  out.push(`  input(非cache)  ${M(st.input).padStart(9)} token  ${pct(st.input, total)}`);
  out.push(`  output          ${M(st.output).padStart(9)} token  ${pct(st.output, total)}`);
  out.push(`  ★合計          ${M(total).padStart(9)} token`);
  out.push('');
  out.push(`  ★1 往復あたりの文脈 ${Math.round(ctx / 1000)}K token`);
  out.push(`  ★★費用 ≒ 文脈 × 往復回数。${Math.round(ctx / 1000)}K × ${st.calls} = ${M(ctx * st.calls)}(実測 ${M(total)})`);
  out.push('');
  const dist = [...st.batchDist.entries()].sort((a, b) => a[0] - b[0]);
  const one = st.batchDist.get(1) || 0;
  out.push(`  tool を呼んだ往復 ${st.reqWithTools} / 呼ばない往復 ${st.noToolCalls}(${pct(st.noToolCalls, st.calls)})`);
  out.push(`  1 往復あたりの tool 本数: ${dist.map(([k, v]) => `${k}本=${v}`).join(' ')}`);
  out.push(`  ★独立呼び出しの束ね率: ${pct(st.reqWithTools - one, st.reqWithTools)}(${st.reqWithTools - one} / ${st.reqWithTools})`);
  return out.join('\n');
}

function toolTable(st) {
  const names = [...new Set([...st.inChars.keys(), ...st.outChars.keys()])];
  const tot = names.reduce((s, n) => s + (st.inChars.get(n) || 0) + (st.outChars.get(n) || 0), 0);
  names.sort((a, b) => ((st.inChars.get(b) || 0) + (st.outChars.get(b) || 0)) - ((st.inChars.get(a) || 0) + (st.outChars.get(a) || 0)));
  const out = ['', '道具別の文脈占有(入力 = assistant が書いた分 / 出力 = 結果。単位 M 字)', ''];
  out.push('道具                          入力M    出力M    合計M      %   呼出   1回入力  1回出力');
  for (const n of names) {
    const i = st.inChars.get(n) || 0, o = st.outChars.get(n) || 0;
    if (i + o < tot / 1000) continue;
    const k = st.inN.get(n) || 0, ko = st.outN.get(n) || 1;
    out.push(`${n.slice(0, 28).padEnd(28)} ${(i / 1e6).toFixed(2).padStart(6)} ${(o / 1e6).toFixed(2).padStart(8)} ${((i + o) / 1e6).toFixed(2).padStart(8)} ${pct(i + o, tot).padStart(7)} ${String(k).padStart(6)} ${String(Math.round(i / (k || 1))).padStart(8)} ${String(Math.round(o / ko)).padStart(8)}`);
  }
  out.push(`合計 ${(tot / 1e6).toFixed(2)}M 字 ≒ ${(tot / 4 / 1e6).toFixed(1)}M token 相当`);
  return out.join('\n');
}

function wasteTable(st) {
  const ctx = st.calls ? st.ctxSum / st.calls : 0;
  const out = ['', '無駄(往復回数を減らせる余地。1 往復 = 文脈 1 回分)', ''];
  out.push(`  ★完全に同じ tool 呼び出しの繰り返し ${st.dup} / ${st.toolCalls}(${pct(st.dup, st.toolCalls)})`);
  const d = [...st.dupBy.entries()].sort((a, b) => b[1] - a[1]).slice(0, 6);
  for (const [n, v] of d) out.push(`      ${n.padEnd(26)} ${v} 回`);
  out.push(`      → ${M(st.dup * ctx)} token 相当`);
  out.push('');
  out.push(`  ★読み取り専用の連続走 ${st.roRuns} 本(読み取り専用 ${st.roCalls} 呼び出し)`);
  out.push(`      束ねれば減らせる往復(上限) ${st.roSaved} → ${M(st.roSaved * ctx)} token 相当`);
  out.push(`      ※上限である。後の呼び出しが前の結果に依存する走も混ざる。`);
  return out.join('\n');
}

function agentTable(st) {
  const rows = [];
  for (const [id, n] of st.notifs) {
    const l = st.launches.get(id);
    const p = l ? l.prompt : '';
    const cjk = (p.match(CJK_RE) || []).length;
    rows.push({
      ts: (n.ts || '').slice(0, 16), tokens: n.tokens, tools: n.tools, ms: n.ms,
      chars: p.length, cjk: p.length ? cjk / p.length : 0, desc: (l ? l.desc : '?').slice(0, 34),
    });
  }
  rows.sort((a, b) => a.ts.localeCompare(b.ts));
  const use = rows.filter((r) => r.tools >= 5);
  const out = ['', `agent ${rows.length} 体(tool_uses >= 5 の ${use.length} 体で統計)`, ''];
  out.push('  ★注: subagent_tokens は cache_read を含まない。「新しく作った量」であって請求量ではない。');
  out.push('');
  out.push('時刻              tokens  tools  brief字  CJK率  desc');
  for (const r of rows) out.push(`${r.ts}  ${String(r.tokens).padStart(7)}  ${String(r.tools).padStart(5)}  ${String(r.chars).padStart(7)}  ${r.cjk.toFixed(2).padStart(5)}  ${r.desc}`);
  if (use.length >= 4) {
    const corr = (f) => {
      const xs = use.map(f), ys = use.map((r) => r.tokens), n = xs.length;
      const mx = xs.reduce((a, b) => a + b, 0) / n, my = ys.reduce((a, b) => a + b, 0) / n;
      let sxy = 0, sxx = 0, syy = 0;
      for (let i = 0; i < n; i++) { const dx = xs[i] - mx, dy = ys[i] - my; sxy += dx * dy; sxx += dx * dx; syy += dy * dy; }
      return sxx && syy ? sxy / Math.sqrt(sxx * syy) : 0;
    };
    out.push('');
    out.push('  subagent_tokens との相関 r:');
    out.push(`    tool_uses     ${corr((r) => r.tools).toFixed(3)}`);
    out.push(`    duration_ms   ${corr((r) => r.ms).toFixed(3)}`);
    out.push(`    brief 字数    ${corr((r) => r.chars).toFixed(3)}`);
    out.push(`    ★CJK 率      ${corr((r) => r.cjk).toFixed(3)}   ← brief の言語`);
    const en = use.filter((r) => r.cjk < 0.02), ja = use.filter((r) => r.cjk >= 0.02);
    if (en.length && ja.length) {
      // tool_uses で調整した残差で比べる
      const xs = use.map((r) => r.tools), ys = use.map((r) => r.tokens), n = xs.length;
      const mx = xs.reduce((a, b) => a + b, 0) / n, my = ys.reduce((a, b) => a + b, 0) / n;
      let sxy = 0, sxx = 0; for (let i = 0; i < n; i++) { sxy += (xs[i] - mx) * (ys[i] - my); sxx += (xs[i] - mx) ** 2; }
      const b1 = sxy / sxx, b0 = my - b1 * mx;
      const res = (r) => r.tokens - (b0 + b1 * r.tools);
      const mean = (a) => a.reduce((x, z) => x + z, 0) / a.length;
      const sd = (a) => (a.length > 1 ? Math.sqrt(a.reduce((s, v) => s + (v - mean(a)) ** 2, 0) / (a.length - 1)) : 0);
      const re = en.map(res), rj = ja.map(res);
      const se = Math.sqrt(sd(re) ** 2 / re.length + sd(rj) ** 2 / rj.length);
      out.push('');
      out.push(`  ★英語 brief(CJK<2%) n=${en.length} と日本語 brief n=${ja.length} を tool_uses で調整して比較:`);
      out.push(`    tokens = ${Math.round(b0)} + ${Math.round(b1)} * tool_uses`);
      out.push(`    残差の差(英−日) ${Math.round(mean(re) - mean(rj))} ± ${Math.round(se)} token  → t = ${se ? ((mean(re) - mean(rj)) / se).toFixed(2) : '—'}`);
    }
  }
  return out.join('\n');
}

// ── selftest ─────────────────────────────────────────────
function selftest() {
  const dir = mkdtempSync(join(tmpdir(), 'sc-'));
  const f = join(dir, 't.jsonl');
  const L = [];
  const asst = (req, blocks, usage) => JSON.stringify({ type: 'assistant', requestId: req, timestamp: '2026-09-06T00:00:00Z', message: { content: blocks, usage } });
  const user = (blocks) => JSON.stringify({ type: 'user', timestamp: '2026-09-06T00:00:00Z', message: { content: blocks } });
  const U = { input_tokens: 1, cache_read_input_tokens: 100, cache_creation_input_tokens: 10, output_tokens: 5 };
  // 1 回の API 呼び出しが 2 レコードに分かれる(thinking + tool_use)。usage は同じ。
  L.push(asst('r1', [{ type: 'thinking', thinking: 'x' }], U));
  L.push(asst('r1', [{ type: 'tool_use', id: 'u1', name: 'Read', input: { file_path: '/a' } }], U));
  L.push(user([{ type: 'tool_result', tool_use_id: 'u1', content: 'AAAA' }]));
  // 同じ Read をもう一度(重複)
  L.push(asst('r2', [{ type: 'tool_use', id: 'u2', name: 'Read', input: { file_path: '/a' } }], U));
  L.push(user([{ type: 'tool_result', tool_use_id: 'u2', content: 'AAAA' }]));
  // 束ねた 2 本(同一 requestId に tool_use 2 つ)
  L.push(asst('r3', [
    { type: 'tool_use', id: 'u3', name: 'Grep', input: { pattern: 'p' } },
    { type: 'tool_use', id: 'u4', name: 'Bash', input: { command: 'echo hi' } },
  ], U));
  // agent の起動と完了
  L.push(JSON.stringify({ type: 'user', timestamp: '2026-09-06T00:01:00Z', toolUseResult: { agentId: 'ag1', description: 'd', prompt: 'hello world' }, message: { content: 'x' } }));
  L.push(JSON.stringify({ type: 'queue-operation', timestamp: '2026-09-06T00:02:00Z', content: '<task-notification><task-id>ag1</task-id><usage><subagent_tokens>1234</subagent_tokens><tool_uses>7</tool_uses><duration_ms>900</duration_ms></usage></task-notification>' }));
  writeFileSync(f, L.join('\n'), 'utf8');
  return scan([f]).then((st) => {
    const t = [];
    const ck = (name, got, want) => t.push({ name, ok: got === want, got, want });
    ck('requestId で重複排除した往復数', st.calls, 3);
    ck('cache_read の合計', st.cacheRead, 300);
    ck('output の合計', st.output, 15);
    ck('tool 呼び出し総数', st.toolCalls, 4);
    ck('完全に同じ呼び出しの重複', st.dup, 1);
    ck('tool を呼んだ往復', st.reqWithTools, 3);
    ck('2 本束ねた往復', st.batchDist.get(2) ?? 0, 1);
    // Read, Read, Grep が連続し(`echo hi` は読み取り専用でない)、走の長さ 3 → 減らせるのは 2
    ck('読み取り専用の連続走で減らせる往復', st.roSaved, 2);
    ck('agent の数', st.notifs.size, 1);
    ck('agent の tool_uses', st.notifs.get('ag1')?.tools, 7);
    ck('brief を突き合わせた', st.launches.get('ag1')?.prompt, 'hello world');
    rmSync(dir, { recursive: true, force: true });
    const ng = t.filter((x) => !x.ok);
    for (const x of t) if (!x.ok) console.log(`  NG ${x.name}: got ${x.got} want ${x.want}`);
    console.log(`selftest ${t.length - ng.length}/${t.length}`);
    return ng.length ? 1 : 0;
  });
}

// ── main ─────────────────────────────────────────────────
if (has('--selftest')) {
  process.exit(await selftest());
}

let files = vals('--file').filter(Boolean);
if (!files.length) {
  const dir = findProjectDir(val('--dir'));
  if (!dir) { console.error('セッション記録の置き場が見つからない。--dir か --file で指定する。'); process.exit(2); }
  files = readdirSync(dir).filter((n) => n.endsWith('.jsonl')).map((n) => join(dir, n));
  if (!files.length) { console.error(`${dir} に .jsonl が無い。`); process.exit(2); }
  if (!has('--json')) console.log(`置き場: ${dir}\n対象: ${files.map((f) => basename(f)).join(', ')}\n`);
}

const st = await scan(files);
if (has('--json')) {
  const total = st.input + st.cacheRead + st.cacheCreate + st.output;
  console.log(JSON.stringify({
    calls: st.calls, total, cacheRead: st.cacheRead, cacheCreate: st.cacheCreate, output: st.output,
    ctxPerCall: st.calls ? Math.round(st.ctxSum / st.calls) : 0,
    toolCalls: st.toolCalls, batched: st.reqWithTools - (st.batchDist.get(1) || 0),
    dup: st.dup, roRuns: st.roRuns, roSaved: st.roSaved,
    agents: [...st.notifs].map(([id, n]) => {
      const p = st.launches.get(id)?.prompt ?? '';
      return { id, tokens: n.tokens, tools: n.tools, ms: n.ms, briefChars: p.length, cjk: p.length ? +((p.match(CJK_RE) || []).length / p.length).toFixed(3) : null };
    }),
  }, null, 1));
} else {
  const all = has('--all');
  console.log(summary(st));
  if (all || has('--tools')) console.log(toolTable(st));
  if (all || has('--waste')) console.log(wasteTable(st));
  if (all || has('--agents')) console.log(agentTable(st));
  if (!all && !has('--tools') && !has('--waste') && !has('--agents')) {
    console.log('\n(--tools 道具別 / --waste 無駄 / --agents agent 別 / --all 全部)');
  }
}
