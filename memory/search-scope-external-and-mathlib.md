---
name: search-scope-external-and-mathlib
description: 既定の grep は external/ を .ignore で拾うが mathlib は拾わない。mathlib は必ず明示パスで検索する。
metadata:
  type: feedback
---

★2026-08-19 実測。`D:\Math_ABC3` の検索範囲は 3 層に分かれる。

| 対象 | 既定の `Grep` / `rg` | 理由 |
|---|---|---|
| `lean/ABC3/`(自前) | **入る** | 追跡されている |
| `external/`(第三者 Lean) | **入る** | `.gitignore` で外し、`.ignore` の `!external/` で打ち消した |
| `lean/.lake/packages/mathlib/` | **入らない** | `lean/.lake/` がディレクトリごと無視されている |

★★**mathlib を再包含することはできない**。gitignore の意味論で
「親ディレクトリが除外されていると中身は再包含できない」ため、
`!lean/.lake/packages/*/Mathlib/` は効かない(実測済み)。
`!lean/.lake/` まで戻すとビルド成果物が数 GB 入るので、やってはいけない。

**How to apply:**

- mathlib を調べるときは**必ず明示パスを渡す**:
  `Grep(pattern, path: "lean/.lake/packages/mathlib/Mathlib/…")`
- `external/` は既定で入るので、「この定理は既に誰かが書いていないか」は
  ルートから普通に grep すればよい。
- 逆に**自前のコードだけを見たいとき**は `path: "lean/ABC3"` を渡す。
  そうしないと `external/` の同名宣言が混ざる。

**Why:** `external/` を落としても検索に入らなければ意味がない。
`.gitignore` だけだと入らず、`.ignore`(ripgrep 専用、git は読まない)で
打ち消して初めて既定検索に入る。★git 側は無視のままなので、
ABC3b の `git add -A` にも巻き込まれない。

関連: [[measure-mathlib-before-skeleton]] [[parallel-session-sweeps-my-files]]
