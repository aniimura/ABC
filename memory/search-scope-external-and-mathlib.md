---
name: search-scope-external-and-mathlib
description: 既定の grep で mathlib まで届く。external/_refs に .lean を複製し .ignore で打ち消してある。ズレ検出は check-external-refs.sh。
metadata:
  type: feedback
---

★2026-08-19 に整備。`D:\Math_ABC3` のルートから普通に grep すると、
**自前コード・第三者 Lean プロジェクト・mathlib のすべてに届く**。

| 対象 | 場所 | 既定の `Grep` / `rg` |
|---|---|---|
| 自前 | `lean/ABC3/` | 入る |
| 第三者プロジェクト | `external/{FLT,formal-conjectures,LeanBridge,iut-lean,lean-poly-abc}` | 入る |
| mathlib / batteries / PrimeNumberTheoremAnd | `external/_refs/` | 入る(**複製**) |
| ビルド成果物・原典 PDF | `lean/.lake/`, `ResearchPaper/0_Source/` | 入らない |

## ★仕組みと、そこに至るまでに潰した選択肢

`.gitignore` で外したものは rg も見ない。そこで **ripgrep 専用の `.ignore`**
(git は読まない)に `!external/` を書いて打ち消してある。

★★mathlib は `lean/.lake/` が**親ごと**除外されているため、
gitignore の意味論で「親が除外されていると中身は再包含できない」。以下は実測で全部駄目だった:

- `.ignore` に `!lean/.lake/packages/*/Mathlib/` —— 効かない
- ジャンクション `mklink /J` —— ripgrep は `--follow` なしで辿らない
- `!lean/.lake/` まで戻す —— ビルド成果物が数 GB 入るので不可

★残ったのは**複製**だけ。`tools/sync-external-refs.sh` が `.lean` だけを
`external/_refs/` へ複製する(8173+175+210 ファイル、111MB)。
`external/` は `.gitignore` のままなので、**git には入らず再配布もしない**。

## ★★陳腐化への備え(ここが要点)

複製は必ず古くなる。`sync-external-refs.sh` は `lake-manifest.json` の rev を
`external/_refs/STAMP.json` に焼く。**`tools/check-external-refs.sh` が rev を比較**し、
ズレていたら NG を返す(数秒)。

**How to apply:**

1. `lake update` をした直後、あるいは mathlib の有無を根拠に判断する前に
   `bash tools/check-external-refs.sh` を走らせる。ズレていたら sync を実行する。
2. **自前のコードだけを見たいときは `path: "lean/ABC3"` を渡す**。
   渡さないと mathlib や FLT の同名宣言が混ざる。
3. 検索コストは気にしなくてよい —— 122MB / 10600 ファイルで **0.1〜0.2 秒**(実測)。

**Why:** `external/` を落としても検索に入らなければ意味がない。
そして「mathlib に無い」という判断こそ最も高くつく間違いなので、
mathlib が既定で入ることの価値が複製のコストを上回る。

関連: [[measure-mathlib-before-skeleton]] [[parallel-session-sweeps-my-files]]
