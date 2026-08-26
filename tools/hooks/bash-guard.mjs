#!/usr/bin/env node
// ABC3 —— Bash の PreToolUse ガード
//
// 2026-08-21 の実測で出た 3 つの損失を潰す:
//   R1  `lake env lean <file>`(数分)  → MCP `lean_check`(0.01 秒)へ誘導して deny
//   R2  作業ディレクトリのドリフト     → `cd /d/Math_ABC3/lean &&` を前置(自動書き換え)
//   R3  mathlib 全体への grep/rg      → deny(`#check` / `#print` で聞く)
//
// 依存なし。標準入力に PreToolUse の JSON が来る。
// ★迷ったときは**通す**(何も出力しない)。ガードが作業を止めないことを優先する。
//
// ★逃げ道は 2 つ:
//   * コマンドに `#no-guard`    …… すべての規則を素通り
//   * コマンドに `#full-check`  …… R1 だけ素通り
//
// ★★誤爆の実例(2026-08-21、導入初日):
//   コミットメッセージのヒアドキュメントに `lake env lean` と書いただけで R1 が発火した。
//   → `stripHeredocs` でヒアドキュメントの中身を判定から外した。

const LEAN_DIR_POSIX = '/d/Math_ABC3/lean';

// ★ヒアドキュメントの中身は「コマンド」ではないので、判定から外す。
function stripHeredocs(src) {
  const lines = src.split('\n');
  const out = [];
  let tag = null;
  for (const line of lines) {
    if (tag !== null) {
      if (line.trim() === tag) tag = null;
      continue;
    }
    out.push(line);
    const m = line.match(
      /<<-?\s*(?:'([A-Za-z_][A-Za-z0-9_]*)'|"([A-Za-z_][A-Za-z0-9_]*)"|([A-Za-z_][A-Za-z0-9_]*))/
    );
    if (m) tag = m[1] || m[2] || m[3];
  }
  return out.join('\n');
}

let raw = '';
process.stdin.setEncoding('utf8');
process.stdin.on('data', (c) => { raw += c; });
process.stdin.on('end', () => {
  let payload;
  try { payload = JSON.parse(raw); } catch { process.exit(0); }
  const cmd = payload && payload.tool_input && payload.tool_input.command;
  if (typeof cmd !== 'string' || cmd.length === 0) process.exit(0);
  if (cmd.indexOf('#no-guard') >= 0) process.exit(0);

  const scan = stripHeredocs(cmd);

  const deny = (reason) => {
    process.stdout.write(JSON.stringify({
      hookSpecificOutput: {
        hookEventName: 'PreToolUse',
        permissionDecision: 'deny',
        permissionDecisionReason: reason,
      },
    }));
    process.exit(0);
  };
  const rewrite = (newCmd, note) => {
    const updated = Object.assign({}, payload.tool_input, { command: newCmd });
    process.stdout.write(JSON.stringify({
      hookSpecificOutput: {
        hookEventName: 'PreToolUse',
        updatedInput: updated,
      },
      systemMessage: note,
    }));
    process.exit(0);
  };

  // ── R1: `lake env lean <file>` ──────────────────────────────────────────
  if (/\blake\s+env\s+lean\b/.test(scan) && scan.indexOf('#full-check') < 0) {
    deny(
      'ABC3 ガード R1: `lake env lean` はファイル全体を再検査するので 1 往復が数分になります。\n'
      + '★代わりに MCP の `mcp__abc3-lean__lean_check`(基準環境の再利用、0.01〜0.03 秒)を使ってください。\n'
      + '  まだ起きていなければ `mcp__abc3-lean__lean_start(["ABC3.Found"])`(初回のみ ~120 秒)。\n'
      + '★ディスクに書いたあとの最終確認は `lake build`(こちらは通ります)。\n'
      + '★どうしても要るときはコマンドに `#full-check` を含めてください。'
    );
  }

  // ── R3: mathlib 全体への grep/rg ───────────────────────────────────────
  if (/\b(grep|rg)\b/.test(scan) && /[.]lake.packages/.test(scan)) {
    deny(
      'ABC3 ガード R3: `.lake/packages`(mathlib 全体)への grep は分単位でかかり、たいてい空振りします。\n'
      + '★補題名・定義を知りたいだけなら `mcp__abc3-lean__lean_check` に\n'
      + '  `#check @Foo` / `open X in #print Foo` を 10 個並べて 1 回で聞いてください(0.01 秒)。'
    );
  }

  // ── R2: 作業ディレクトリの前置 ──────────────────────────────────────────
  const startsWithCd = /^\s*cd\s/.test(cmd);
  // ★リポジトリ直下のパスを触るコマンドは、cwd を動かすと壊れる(2026-08-21 に誤爆)。
  const usesRepoRoot = /(^|[\s;&|("])(tools|ResearchPaper|lean)[/]/.test(scan);
  const usesLake = /(^|[\s;&|(])lake\s/.test(scan);
  const usesRelAbc3 = /(^|[\s;&|("])ABC3[/]/.test(scan);
  if (!startsWithCd && !usesRepoRoot && (usesLake || usesRelAbc3)) {
    rewrite(
      'cd ' + LEAN_DIR_POSIX + ' && ' + cmd,
      'ABC3 ガード R2: 作業ディレクトリを ' + LEAN_DIR_POSIX + ' に固定しました(相対パスの取り違えを防ぐため)。'
    );
  }

  process.exit(0);
});
