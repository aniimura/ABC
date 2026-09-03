---
name: sibling-session-coordination-via-listagents
description: 共有ワークツリーの並行セッション(ABC3b・math-abc3-d3等)には ListAgents/SendMessage で直接連絡できる。着手前にトラックを聞くと衝突を防げる。
metadata:
  type: feedback
---

★2026-09-03、ABC3c 起動時に実測。`D:\Math_ABC3` を共有する他の Claude セッションは
**同一マシン上の別セッションとして `ListAgents` に列挙される**(名前は各セッションが
`/rename` で付けた表示名。例: `ABC3b`・`math-abc3-d3`)。`SendMessage({to: "<名前>", ...})`
で直接メッセージを送れ、相手が手すきになったタイミングで返信が
`<cross-session-message from="...">` として届く。

**やったこと**: 起動直後に `ABC3b` と `math-abc3-d3` へ「今どのファイル/トラックを
触っているか」を尋ねた。`math-abc3-d3` から即座に返信があり、
「`lean/ABC3/**` は触っていない、`dependency-graph.html` 生成ツールと
`papers.json` 登記をしている最中」と分かった。`ResearchPaper/mathlib-gap.json` を見ると
`ABC3b` は GenEll `Lemma 3.5` に数百ノード単位で数日がかりで集中していることも分かった
(git のコミット不要、ファイルを読むだけで検知できた)。

**Why:** [[parallel-session-sweeps-my-files]] は「衝突したときの後始末」だが、
これは**衝突を先に避ける**手段。CLAUDE.md の未コミット差分(他セッションが編集中)を見て
「今まさに誰かが触っている」と気づけたのも大きい——`git status`/`git diff` で
未コミットの変更を見れば、共有ワークツリー上の**今の**activityが分かる。

**How to apply:**
- 新しいセッションを起動したら、着手前に必ず `ListAgents` で同名プロジェクトの
  他セッション(busy/idle)を確認する。
- busy な相手には一言メッセージで「今のトラック」を聞く。返信を待つ間、
  `git status --short` / 各 `Found/<Track>/` の最新 mtime / `ResearchPaper/*.json` の
  最近のキーで**間接的に**現在地を推定できる(相手の返信が無くても手が止まらない)。
- GenEll は [[genell-track-b]] のとおり ABC3b の担当。衝突を避けるなら
  `lean/ABC3/**/GenEll/**` と `ResearchPaper/genell-goal.md` 系には触らない。
- 相手の返信で「良い」と言われたファイル(例: `papers.json`)以外は
  自由に着手してよい——ただし [[parallel-session-sweeps-my-files]] の
  コミット規律(1 ファイル通ったら即コミット、`git add -A` 禁止)は変わらず適用する。

関連: [[parallel-session-sweeps-my-files]] / [[child-session-usage]] / [[genell-track-b]] /
[[arakelov-galois-skeleton-counts]] / [[leaf-first-with-graph-feedback]]
