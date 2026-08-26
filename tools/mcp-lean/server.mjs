#!/usr/bin/env node
// ABC3 Lean REPL MCP —— 宣言 1 個を、ビルド済み import に対して検査する。
//
// ★動機(2026-08-17 実測):
//   `lake env lean ABC3/Found/FrdI/Prop32Frob.lean`(3,804 行) …… 数分 / 1 回
//   `lake env lean ABC3/Found/FrdI/Example43.lean`(322 行)   …… 9 秒 / 1 回
//   本サーバ(env 再利用)                                      …… 0.01〜0.03 秒 / 1 回
//   ★import の読み込みは**セッション中 1 回だけ**(約 90 秒)。
//
// 依存なし(Node 22 の標準ライブラリのみ)。MCP は stdio + 行区切り JSON-RPC 2.0。

import { spawn, spawnSync } from 'node:child_process';
import { createInterface } from 'node:readline';
import { fileURLToPath } from 'node:url';
import path from 'node:path';

// このファイルは <repo>/tools/mcp-lean/server.mjs にある。
const HERE = path.dirname(fileURLToPath(import.meta.url));
const REPO = path.resolve(HERE, '..', '..');

const LEAN_DIR = process.env.ABC3_LEAN_DIR ?? path.join(REPO, 'lean');
const REPL_BIN =
  process.env.ABC3_LEAN_REPL ??
  path.join(REPO, 'tools', 'lean-repl', '.lake', 'build', 'bin',
    process.platform === 'win32' ? 'repl.exe' : 'repl');
const DEFAULT_TIMEOUT_MS = Number(process.env.ABC3_LEAN_TIMEOUT_MS ?? 600000);

// ★2026-08-20 の恒久対処: Lean REPL は **環境スナップショットを一切解放しない**。
// `lean_check` を 1 回呼ぶたびに新しい `Environment` が REPL 内部の配列に積まれ、
// 長いセッションでは 18〜46 GB まで膨らむ(実測)。
// そこで **一定回数 / 一定メモリでプロセスを建て直す**。
// `addToEnv` で積んだ宣言は `state.envLog` に控えてあるので、建て直しのあとに
// 再生する——呼び出し側から見て基準環境は失われない。
// ★★回数の方が**主たる歯止め**である。`tasklist` が返すのは working set なので、
// メモリ圧が掛かって OS にトリムされると 17 GB のプロセスが 1.7 GB に見える
// (2026-08-20 に実測)。メモリ判定は「常駐して大きいとき」にしか効かないので、
// 決定的な回数の方を主に据える。
const MAX_MB = Number(process.env.ABC3_LEAN_MAX_MB ?? 6000);
const MAX_CHECKS = Number(process.env.ABC3_LEAN_MAX_CHECKS ?? 120);
const MEM_CHECK_EVERY = Number(process.env.ABC3_LEAN_MEM_CHECK_EVERY ?? 20);
// ★★2026-08-20 の追加対処: **遊んでいる REPL を回収する**。
// 上の建て直しは `lean_check` の中でしか走らないので、
// 途中から `lake build` に切り替えたセッションでは
// **使われない REPL が数 GB を抱えたまま永遠に残る**（2026-08-20 に実測）。
// 一定時間触られていなければ落とす——次の `lean_check` で自動的に読み直される。
const IDLE_KILL_MS = Number(process.env.ABC3_LEAN_IDLE_KILL_SEC ?? 900) * 1000;

/** @type {{proc: import('node:child_process').ChildProcess|null, baseEnv: number|null, imports: string[], leanPath: string|null, busy: boolean, buf: string[], pending: ((s: string) => void)|null}} */
const state = {
  proc: null,
  baseEnv: null,
  imports: [],
  leanPath: null,
  busy: false,
  buf: [],
  pending: null,
  checks: 0,
  /** ★`addToEnv` で基準環境に積んだコード片。建て直しのときに再生する。 */
  envLog: [],
  recycles: 0,
  /** ★最後に REPL を触った時刻（遊び回収の判定用）。 */
  lastUse: Date.now(),
};

/** ★生きている `repl.exe` の物理メモリを報告する(2026-08-20 の 46 GB 事故の再発検知)。 */
function replMemory() {
  if (process.platform !== 'win32') return 'repl メモリ: (未計測)';
  try {
    const r = spawnSync(
      'tasklist',
      ['/FI', 'IMAGENAME eq repl.exe', '/FO', 'CSV', '/NH'],
      { encoding: 'utf8', windowsHide: true },
    );
    const rows = (r.stdout ?? '')
      .split(/\r?\n/)
      .map((l) => l.match(/^"repl\.exe","(\d+)",[^,]*,[^,]*,"([\d,. ]+) K"$/))
      .filter(Boolean)
      .map((m) => ({ pid: m[1], mb: Math.round(Number(m[2].replace(/[^\d]/g, '')) / 1024) }));
    if (rows.length === 0) return 'repl メモリ: (repl.exe は動いていない)';
    const total = rows.reduce((a, b) => a + b.mb, 0);
    const detail = rows.map((x) => `pid ${x.pid}: ${x.mb} MB`).join(' / ');
    return `repl メモリ: 合計 ${total} MB (${detail})` +
      (total > MAX_MB ? '  ★★閾値超——次の lean_check で自動的に建て直す' : '') +
      '  ※working set なのでトリムされると小さく見える';
  } catch {
    return 'repl メモリ: (計測に失敗)';
  }
}

/** ★生きている `repl.exe` の物理メモリ合計(MB)。数えられなければ 0。
 * 他セッションの `repl.exe` も含む上界である——安全側に倒すためこれでよい。 */
function replMemoryMB() {
  if (process.platform !== 'win32') return 0;
  try {
    const r = spawnSync(
      'tasklist',
      ['/FI', 'IMAGENAME eq repl.exe', '/FO', 'CSV', '/NH'],
      { encoding: 'utf8', windowsHide: true },
    );
    return (r.stdout ?? '')
      .split(/\r?\n/)
      .map((l) => l.match(/^"repl\.exe","(\d+)",[^,]*,[^,]*,"([\d,. ]+) K"$/))
      .filter(Boolean)
      .reduce((a, m) => a + Number(m[2].replace(/[^\d]/g, '')) / 1024, 0);
  } catch {
    return 0;
  }
}

function log(msg) {
  process.stderr.write(`[mcp-lean] ${msg}\n`);
}

// ★2026-08-20 に踏んだ罠: Windows では `shell: true` で起こすので
// `cmd.exe → lake.exe → lake.exe → repl.exe` の 4 段になる。
// `proc.kill()` は先頭の `cmd.exe` しか殺さないので、**mathlib を抱えた
// `repl.exe`(17〜22 GB)が孤児として残り続ける**。実測でコミットが
// 104 GB / 119 GB まで来ていた。プロセスツリーごと殺すこと。
function killTree(proc) {
  if (!proc || proc.pid == null) return;
  if (process.platform === 'win32') {
    try {
      spawnSync('taskkill', ['/pid', String(proc.pid), '/T', '/F'], {
        stdio: 'ignore',
        windowsHide: true,
      });
    } catch {
      /* ignore */
    }
  }
  try {
    proc.kill();
  } catch {
    /* already gone */
  }
}

function killRepl() {
  killTree(state.proc);
  state.proc = null;
  state.baseEnv = null;
  state.pending = null;
  state.buf = [];
  state.busy = false;
}

/** ★基準環境ごと捨てる(`lean_start` / `lean_reset` 用)。 */
function killReplHard() {
  killRepl();
  state.envLog = [];
  state.checks = 0;
}

function startRepl() {
  killRepl();
  // ★`lake env` 経由で起動する。Windows では `libleanshared.dll` の解決に
  // ツールチェインの PATH が要るので、LEAN_PATH だけ渡す直起動では動かない
  // (2026-08-17 に実測: 直起動は無応答、`lake env` 経由は 90 秒で import 完了)。
  const proc = spawn('lake', ['env', REPL_BIN], {
    cwd: LEAN_DIR,
    shell: process.platform === 'win32',
  });
  proc.on('error', (e) => log(`repl spawn error: ${e.message}`));
  proc.on('exit', (code) => {
    log(`repl exited (code=${code})`);
    if (state.proc === proc) {
      state.proc = null;
      state.baseEnv = null;
    }
  });
  const rl = createInterface({ input: proc.stdout });
  rl.on('line', (line) => {
    if (line.trim() === '' && state.buf.length > 0) {
      const text = state.buf.join('\n');
      state.buf = [];
      const done = state.pending;
      state.pending = null;
      if (done) done(text);
    } else if (line.trim() !== '') {
      state.buf.push(line);
    }
  });
  proc.stderr.on('data', (d) => log(`repl stderr: ${String(d).trim()}`));
  state.proc = proc;
  return proc;
}

/** REPL に 1 コマンド送って応答(JSON 文字列)を待つ。 */
function replSend(obj, timeoutMs) {
  return new Promise((resolve, reject) => {
    if (!state.proc) return reject(new Error('REPL が起動していない'));
    if (state.busy) return reject(new Error('REPL は処理中(直列にしか使えない)'));
    state.busy = true;
    state.lastUse = Date.now();
    const timer = setTimeout(() => {
      state.busy = false;
      killRepl();
      reject(
        new Error(
          `${Math.round(timeoutMs / 1000)} 秒で応答が無いので REPL を落とした。` +
            `次の呼び出しで自動的に再起動して import を読み直す。`,
        ),
      );
    }, timeoutMs);
    state.pending = (text) => {
      clearTimeout(timer);
      state.busy = false;
      resolve(text);
    };
    state.proc.stdin.write(JSON.stringify(obj) + '\n\n');
  });
}

/** imports を読み込んで baseEnv を作る。 */
async function loadImports(imports, timeoutMs) {
  if (!state.proc) startRepl();
  const cmd = imports.map((m) => `import ${m}`).join('\n');
  const t0 = Date.now();
  const text = await replSend({ cmd, env: null }, timeoutMs);
  let res;
  try {
    res = JSON.parse(text);
  } catch {
    throw new Error(`REPL の応答が JSON でない:\n${text.slice(0, 2000)}`);
  }
  const errs = (res.messages ?? []).filter((m) => m.severity === 'error');
  if (errs.length > 0) {
    throw new Error(
      `import に失敗:\n` + errs.map((m) => m.data).join('\n').slice(0, 2000),
    );
  }
  state.baseEnv = res.env ?? null;
  state.imports = imports;
  return { seconds: (Date.now() - t0) / 1000, env: state.baseEnv };
}

/** ★★★**建て直しが要るか**——回数かメモリのどちらかが閾値を超えたら。 */
function needsRecycle() {
  if (state.checks >= MAX_CHECKS) return `検査 ${state.checks} 回`;
  if (state.checks > 0 && state.checks % MEM_CHECK_EVERY === 0) {
    const mb = Math.round(replMemoryMB());
    if (mb >= MAX_MB) return `repl.exe 合計 ${mb} MB`;
  }
  return null;
}

/** ★★★★**REPL を建て直し、`addToEnv` で積んだ宣言を再生する**。
 * 呼び出し側から見て基準環境は保たれる(env の id だけが変わる)。 */
async function recycleRepl(timeoutMs, why) {
  const saved = state.envLog.slice();
  log(`recycle: ${why} —— REPL を建て直す(再生 ${saved.length} 件)`);
  killRepl();
  state.checks = 0;
  state.recycles += 1;
  await loadImports(state.imports, timeoutMs);
  for (const code of saved) {
    const text = await replSend({ cmd: code, env: state.baseEnv }, timeoutMs);
    let res;
    try {
      res = JSON.parse(text);
    } catch {
      throw new Error(`建て直しの再生で REPL の応答が JSON でない:\n${text.slice(0, 1000)}`);
    }
    const errs = (res.messages ?? []).filter((m) => m.severity === 'error');
    if (errs.length > 0 || res.env == null) {
      throw new Error(
        `建て直しの再生に失敗した。lean_start をやり直すこと:\n` +
          errs.map((m) => m.data).join('\n').slice(0, 1000),
      );
    }
    state.baseEnv = res.env;
  }
  state.envLog = saved;
}

function fmtMessages(res) {
  const msgs = res.messages ?? [];
  if (msgs.length === 0) return null;
  return msgs
    .map((m) => {
      const p = m.pos ? `${m.pos.line}:${m.pos.column}` : '?';
      return `${m.severity} ${p}\n${m.data}`;
    })
    .join('\n\n');
}

// ---------------------------------------------------------------- tools

const TOOLS = [
  {
    name: 'lean_start',
    description:
      'Lean REPL を起動し、指定した import を読み込んで基準環境を作る。' +
      '初回のみ 90 秒前後かかる(mathlib + プロジェクトの olean を読む)。' +
      'olean を再ビルドしたあとは、これを呼び直して読み込み直すこと。',
    inputSchema: {
      type: 'object',
      properties: {
        imports: {
          type: 'array',
          items: { type: 'string' },
          description:
            'import するモジュール名の配列。例: ["ABC3.Found.FrdI.Prop32Frob"]',
        },
        timeoutSeconds: { type: 'number' },
      },
      required: ['imports'],
    },
  },
  {
    name: 'lean_check',
    description:
      '基準環境に対して Lean のコード片(定理 1 個など)を検査し、' +
      'エラー・警告・#check の出力を返す。★基準環境は汚さない(独立検査)。' +
      'lean_start がまだなら、直前の imports で自動起動する。' +
      '0.01〜0.03 秒で返るので、ファイル全体の再検査の代わりに使う。',
    inputSchema: {
      type: 'object',
      properties: {
        code: { type: 'string', description: '検査する Lean のコード片' },
        addToEnv: {
          type: 'boolean',
          description:
            'true なら、エラーが無かったときに基準環境へ積む(ファイルを少しずつ組み立てるとき)。既定 false。',
        },
        timeoutSeconds: { type: 'number' },
      },
      required: ['code'],
    },
  },
  {
    name: 'lean_status',
    description: 'REPL の状態(起動しているか・読み込んだ import・基準環境 id)を返す。',
    inputSchema: { type: 'object', properties: {} },
  },
  {
    name: 'lean_reset',
    description:
      'REPL を落として基準環境を捨てる。次の lean_start / lean_check で再起動する。',
    inputSchema: { type: 'object', properties: {} },
  },
];

async function callTool(name, args) {
  const timeoutMs = args?.timeoutSeconds
    ? args.timeoutSeconds * 1000
    : DEFAULT_TIMEOUT_MS;

  if (name === 'lean_start') {
    // ★2026-08-20: 以前はプロセスが生きていれば使い回していたが、
    // **REPL は環境スナップショットを一切解放しない**ので、`import` を
    // 呼ぶたびに mathlib 1 組分がプロセス内に積まれる(実測 46 GB)。
    // 建て直す。温まっていれば 4〜5 秒で戻るので、使い回す利点は無い。
    killReplHard();
    const r = await loadImports(args.imports, timeoutMs);
    state.checks = 0;
    return (
      `起動して import を読み込んだ (${r.seconds.toFixed(1)} 秒)。\n` +
      `imports: ${args.imports.join(', ')}\n基準環境 env=${r.env}`
    );
  }

  if (name === 'lean_check') {
    if (state.baseEnv === null) {
      if (state.imports.length === 0) {
        return 'まだ lean_start を呼んでいない。imports を指定して lean_start を呼ぶこと。';
      }
      await loadImports(state.imports, timeoutMs);
    }
    // ★恒久対処: 積み上がる前に建て直す(`addToEnv` の宣言は再生される)。
    const why = needsRecycle();
    let recycled = '';
    if (why) {
      await recycleRepl(timeoutMs, why);
      recycled = ` / ★REPL を建て直した(${why}、再生 ${state.envLog.length} 件)`;
    }
    const t0 = Date.now();
    const text = await replSend({ cmd: args.code, env: state.baseEnv }, timeoutMs);
    const dt = ((Date.now() - t0) / 1000).toFixed(2);
    state.checks += 1;
    let res;
    try {
      res = JSON.parse(text);
    } catch {
      return `REPL の応答が JSON でない (${dt} 秒):\n${text.slice(0, 4000)}`;
    }
    const msgs = fmtMessages(res);
    const errs = (res.messages ?? []).filter((m) => m.severity === 'error');
    if (args.addToEnv && errs.length === 0 && res.env != null) {
      state.baseEnv = res.env;
      state.envLog.push(args.code);
    }
    const head =
      errs.length === 0
        ? `OK (${dt} 秒)${args.addToEnv ? ` / 基準環境を env=${state.baseEnv} に更新` : ''}${recycled}`
        : `エラー ${errs.length} 件 (${dt} 秒)${recycled}`;
    const sorries = (res.sorries ?? []).length;
    const tail = sorries > 0 ? `\n\n★sorry ${sorries} 件` : '';
    return msgs ? `${head}\n\n${msgs}${tail}` : `${head}${tail}`;
  }

  if (name === 'lean_status') {
    return (
      `起動: ${state.proc ? 'あり' : 'なし'}\n` +
      `imports: ${state.imports.join(', ') || '(なし)'}\n` +
      `基準環境: ${state.baseEnv ?? '(なし)'}\n` +
      `lean_check 回数: ${state.checks} / ${MAX_CHECKS}(超えたら自動で建て直す)\n` +
      `再生用に控えた宣言: ${state.envLog.length} 件 / 建て直し ${state.recycles} 回\n` +
      `${replMemory()}\n` +
      `自動建て直しの閾値: ${MAX_MB} MB(${MEM_CHECK_EVERY} 回ごとに計測)\n` +
      `LEAN_DIR: ${LEAN_DIR}\nREPL: ${REPL_BIN}`
    );
  }

  if (name === 'lean_reset') {
    killReplHard();
    return 'REPL を落とした(再生用の控えも捨てた)。';
  }

  throw new Error(`未知のツール: ${name}`);
}

// ---------------------------------------------------------------- MCP

function send(msg) {
  process.stdout.write(JSON.stringify(msg) + '\n');
}

const rl = createInterface({ input: process.stdin });
rl.on('line', async (line) => {
  const s = line.trim();
  if (s === '') return;
  let req;
  try {
    req = JSON.parse(s);
  } catch {
    return;
  }
  const { id, method, params } = req;

  try {
    if (method === 'initialize') {
      send({
        jsonrpc: '2.0',
        id,
        result: {
          protocolVersion: params?.protocolVersion ?? '2024-11-05',
          capabilities: { tools: {} },
          serverInfo: { name: 'abc3-lean', version: '0.1.0' },
        },
      });
      return;
    }
    if (method === 'notifications/initialized' || method === 'initialized') return;
    if (method === 'tools/list') {
      send({ jsonrpc: '2.0', id, result: { tools: TOOLS } });
      return;
    }
    if (method === 'tools/call') {
      const text = await callTool(params.name, params.arguments ?? {});
      send({
        jsonrpc: '2.0',
        id,
        result: { content: [{ type: 'text', text: String(text) }] },
      });
      return;
    }
    if (method === 'ping') {
      send({ jsonrpc: '2.0', id, result: {} });
      return;
    }
    if (id !== undefined) {
      send({
        jsonrpc: '2.0',
        id,
        error: { code: -32601, message: `未対応の method: ${method}` },
      });
    }
  } catch (e) {
    if (id !== undefined) {
      send({
        jsonrpc: '2.0',
        id,
        result: { content: [{ type: 'text', text: `エラー: ${e.message}` }], isError: true },
      });
    }
  }
});

// ★2026-08-20 の恒久対処: どの経路で終わっても `repl.exe` を孤児にしない。
// MCP サーバは stdin が閉じたら終わる——これを拾わないと、セッションを閉じても
// `lake.exe → repl.exe` が数十 GB を抱えたまま残る。
function shutdown(code) {
  killRepl();
  process.exit(code ?? 0);
}
// ★★★**遊んでいる REPL の回収**——`lean_check` を呼ばないままになった
// セッションでもメモリを抱えっぱなしにしない。基準環境は `state.envLog` に
// 控えてあるので、次の `lean_check` が import を読み直して再生する。
const idleTimer = setInterval(() => {
  if (!state.proc || state.busy) return;
  if (Date.now() - state.lastUse < IDLE_KILL_MS) return;
  const mb = Math.round(replMemoryMB());
  log(`idle: ${Math.round((Date.now() - state.lastUse) / 1000)} 秒触られていないので REPL を落とす (${mb} MB)`);
  killRepl();
}, 60_000);
idleTimer.unref?.();

process.on('exit', killRepl);
process.on('SIGINT', () => shutdown(0));
process.on('SIGTERM', () => shutdown(0));
process.on('SIGHUP', () => shutdown(0));
process.on('uncaughtException', (e) => {
  log(`uncaught: ${e?.message ?? e}`);
  shutdown(1);
});
rl.on('close', () => shutdown(0));
process.stdin.on('close', () => shutdown(0));
