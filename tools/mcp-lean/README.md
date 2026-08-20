# abc3-lean —— Lean を「宣言 1 個ずつ」検査する MCP

## なぜ作ったか（2026-08-17 の実測）

形式化の作業は「書く → 検査 → エラーを読む → 直す」の往復で、**律速は検査**である。
`lake env lean` は**ファイル全体を再検査する**ので、ファイルが育つほど 1 往復が重くなる。

| 検査のしかた | 1 往復 |
|---|---|
| `lake env lean ABC3/Found/FrdI/Prop32Frob.lean`（3,804 行） | **数分** |
| `lake env lean ABC3/Found/FrdI/Example43.lean`（322 行） | **9 秒** |
| 本 MCP の `lean_check`（基準環境を再利用） | **0.01〜0.02 秒** |

`import` の読み込みだけはセッション中 1 回必要で、**冷たいときで約 90 秒・温まっていれば約 8 秒**。
それ以降の検査は定数時間になる。

★ファイルを小さく割るのも同じ方向の対策だが、本 MCP は**割らなくても**同じ効果を出す
（しかも 9 秒 → 0.01 秒までさらに縮む）。

## 中身

- `tools/lean-repl/` …… [leanprover-community/repl](https://github.com/leanprover-community/repl) の
  `v4.31.0`（本プロジェクトの `lean-toolchain` と一致）。`git` からは無視している。
- `tools/mcp-lean/server.mjs` …… 依存なしの MCP サーバ（Node 22 の標準ライブラリのみ）。

## 使い方

`.mcp.json` に登録済み。**セッションを開き直すと使えるようになる。**

| ツール | 用途 |
|---|---|
| `lean_start(imports)` | REPL を起こして `import` を読む。基準環境ができる |
| `lean_check(code, addToEnv?)` | 基準環境に対してコード片を検査する。既定では基準環境を汚さない |
| `lean_status()` | 起動状態・読み込んだ `import`・基準環境 id |
| `lean_reset()` | 落として捨てる |

典型的な流れ:

1. `lean_start(["ABC3.Found.FrdI.Prop32"])`
2. 定理を 1 個ずつ `lean_check` で通す（0.01 秒）
3. 通ったものをファイルへ書き、最後に一度だけ `lake build` で確定する

★**`lake build` を省いてはならない。** REPL は `olean` を読むだけなので、
ファイルに書いた順序・`variable` の効き方・リンタは `lake build` でしか確かめられない。

## メモリ（2026-08-20 に測った）

`repl.exe` は **mathlib と ABC3 の `olean` を常駐**させるので、起動直後で数 GB。
そのうえ `lean_check(addToEnv: true)` は**環境スナップショットを 1 個ずつ積む**。
プロセスが生きているあいだ、積んだ環境は解放されない。

★実測: 1 セッションで env を 90 個ほど積んだ REPL が **14 GB** に達した。
重いのは個々の証明ではなく、**長生きした REPL の積み上げ**である。

### ★★★★★★恒久対処——自動リサイクル（2026-08-20 に入れた）

**運用の注意で抑えるのはやめた。**サーバが自分で建て直す。

| 仕掛け | 既定値 | 環境変数 |
|---|---|---|
| `lean_check` の回数で建て直す | **120 回** | `ABC3_LEAN_MAX_CHECKS` |
| `repl.exe` の合計メモリで建て直す | **6000 MB** | `ABC3_LEAN_MAX_MB` |
| メモリを測る間隔 | 20 回ごと | `ABC3_LEAN_MEM_CHECK_EVERY` |

★★**`addToEnv` で積んだ宣言は失われない。** サーバは積んだコード片を `state.envLog`
に控えていて、建て直しのあとに**順に再生する**。呼び出し側から見えるのは
`env` の id が変わることと、その 1 回だけ数秒余分に掛かることだけである。

実測（`ABC3_LEAN_MAX_CHECKS=2` で確認）:

```
--- 3 (ここで建て直しが起きるはず) ---
[mcp-lean] recycle: 検査 2 回 —— REPL を建て直す(再生 2 件)
OK (0.00 秒) / 基準環境を env=3 に更新 / ★REPL を建て直した(検査 2 回、再生 2 件)
--- 4 (再生できていれば foo3 が見える) ---
OK (0.02 秒)   info 1:0   3
```

★**回数の方が主たる歯止め**である。`tasklist` が返すのは working set なので、
メモリ圧が掛かって OS にトリムされると **17 GB のプロセスが 1.7 GB に見える**
（2026-08-20 に実測）。メモリ判定は「常駐して大きいとき」にしか効かない。

### ★★★★終了時にツリーごと落とす

セッションを閉じたときに `repl.exe` を孤児にしないため、
`exit` / `SIGINT` / `SIGTERM` / `SIGHUP` / `uncaughtException` /
**stdin の close**（MCP サーバはこれで終わる）のすべてで `killTree` を呼ぶ。

★実測: これを入れる前は、親の `claude.exe` が生きていても `repl.exe` が
17 GB を抱えたまま残った。入れたあとはドライバ終了後に `repl.exe` が 0 個になる。

### ★それでも残る運用の勘所

* 使い捨ての実験は **`addToEnv` を付けない**（既定で基準環境を汚さない）。
  後続が使う補題だけ `addToEnv: true` にする——再生の対象が減って建て直しが速い。
* `lean_status` が回数・控えた宣言の数・メモリ・建て直し回数を報告する。

## 作り直すとき

```
cd tools/lean-repl && lake build      # repl.exe を作る
node --check tools/mcp-lean/server.mjs
```

## 実装上の注意（踏んだ罠）

★Windows では `repl.exe` を**直に spawn してはならない**。`libleanshared.dll` の解決に
ツールチェインの `PATH` が要る。直起動は無応答のまま 600 秒でタイムアウトした。
`lake env <repl>` 経由で起動すること（`server.mjs` にコメントで残してある）。

### ★★★メモリが 46 GB まで膨らんだ（2026-08-20 に踏んだ）

ユーザに「18 GB ほど消費しているが想定通りか」と指摘されて発覚した。
実測は `repl.exe` 1 本で**作業セット 17 GB・コミット 46 GB**、
システム全体でコミット **104.8 GB / 119.7 GB**、空き物理メモリ 7.3 GB、
Memory Compression が 14 GB。原因は**独立した 2 つのバグ**であった。

| # | バグ | 直し方 |
|---|---|---|
| 1 | `killRepl` が**孫を殺していなかった** | `taskkill /T /F` でツリーごと殺す |
| 2 | `lean_start` が**プロセスを使い回していた** | 毎回建て直す |

**(1)** Windows では `shell: true` で起こすので
`cmd.exe → lake.exe → lake.exe → repl.exe` の 4 段になる。
`proc.kill()` は先頭の `cmd.exe` しか殺さないので、
**mathlib を抱えた `repl.exe` が孤児として残り続ける**。
`lean_reset` を呼んでも 17 GB のプロセスが生き残っていた。

**(2)** `loadImports` は `if (!state.proc) startRepl()` だったので、
プロセスが生きていれば**同じプロセスに mathlib をもう一組読み込む**。
★★**REPL は環境スナップショット（`env=0,1,2,…`）を一切解放しない**ので、
`import` のたびに mathlib 1 組分がプロセス内に積まれる。
env=151 まで積んだ結果が 46 GB であった。
建て直しは温まっていれば 4〜5 秒で戻るので、使い回す利点は無い。

★孤児を落としただけで **48 GB 解放された**
（空き 7.3 GB → 39.9 GB、コミット 104.8 GB → 56.7 GB）。

★★再発検知のため `lean_status` が `repl.exe` の物理メモリと
`lean_check` の回数を報告するようにした。

★★★**2026-08-20（同日、2 度目の再発を受けて）恒久対処を入れた**——
上の「自動リサイクル」を見ること。ユーザから 3 度（18 GB / 13.6 GB / 17 GB）
指摘を受けたので、運用の注意で抑えるのをやめてサーバ側で自動化した。

★★★**孤児が出たときの手当て**（`server.mjs` を直す前のセッションでも使える）:

```powershell
# 素性を確かめる（並行セッションのものを巻き込まないこと）
Get-CimInstance Win32_Process -Filter "Name='repl.exe' OR Name='lake.exe'" |
  Select-Object ProcessId, ParentProcessId, Name, CreationDate, CommandLine | Format-List
# 親（node tools/mcp-lean/server.mjs）が死んでいるものが孤児
Stop-Process -Id <repl の pid>,<lake の pid 2 つ> -Force
```
