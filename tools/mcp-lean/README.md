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

## 作り直すとき

```
cd tools/lean-repl && lake build      # repl.exe を作る
node --check tools/mcp-lean/server.mjs
```

## 実装上の注意（踏んだ罠）

★Windows では `repl.exe` を**直に spawn してはならない**。`libleanshared.dll` の解決に
ツールチェインの `PATH` が要る。直起動は無応答のまま 600 秒でタイムアウトした。
`lake env <repl>` 経由で起動すること（`server.mjs` にコメントで残してある）。
