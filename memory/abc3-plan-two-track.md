---
name: abc3-plan-two-track
description: Math_ABC3 の北極星(IUT全体+abc)と、2トラック構成・statementゲート・飛躍の扱いという計画の骨子。
metadata: 
  node_type: memory
  type: project
  originSessionId: fbb27ed5-3b09-4626-a97a-bdf85c79f517
  modified: 2026-08-14T01:53:27.388Z
---

北極星: **望月新一氏のIUT全体(abc予想の証明を含む)のLean形式化**。原典に飛躍がありうることを前提に置き、その局面の姿勢は `D:\Math_ABC3\idea2.md`(Wiles流「目的には頑固、方法には柔軟」)に従う。計画本体は `D:\Math_ABC3\PLAN.md`。

計画の骨子(2026-08-14に検討・実測に基づく):

- **2トラック**: 上から `Skeleton/`(statementのみ、`sorry`)、下から `Found/`(実装、sorry無し)、接続点が `Interface/`(未構築の基礎を `structure` で受ける)。
- **`axiom` 禁止**: `axiom` は検出不能な破滅、`structure` は検出可能な空虚。Lean で実演確認済み。
- **statementゲート G1-G5**: G1出典(PDFページ+逐語+記法)/ G2非空虚witness / G3負の対照 / G4 axiomゼロ / G5弱化禁止。load-bearing set にのみ全部課し、他は G1・G4 のみ。
- **飛躍は追加仮説として型に出す**(`Gap/`)。「証明できない」ではなく `theorem p (h : Data) (H : Gap h) : C`。①モデル化の誤り ②未構築の数学 ③原典側の飛躍 を区別し、**既定は①**、③には falsifier 必須。

**Why:** Lean は `theorem foo : P := by sorry` について型検査しかしない。実演で確認: 条件を1つ取り違えた `structure` は skeleton がビルドを通るだけでなく、`sorry` が「埋まって」しまい、sorry無し・**公理依存ゼロ**で `0 = 1` まで証明できる。skeleton方式は検査されない側(statement)に内容の100%を置くので、statement専用ゲートが要る。

**How to apply:** 新しい skeleton を書くときは必ず G1・G4 を満たす。`Interface` を置いたことを「形式化した」と呼ばない。詳細と再現手順は `PLAN.md` §1。前提の扱いは [[abc3-independent-of-prior-generations]]、争点の所在は [[abc3-disputed-locus-cor312]]。
