---
name: no-wall-decompose-instead
description: 工数の山を「壁」と呼ばない。到達不能と報告する前に、塊を小目標の DAG に割って葉を出す。割ると在庫の測り違いが出てくる。
metadata:
  type: feedback
---

★利用者が 2026-08-18 に CLAUDE.md へ追記した行:

> 姿勢 : 工数の山を「壁」と呼ばない。既知数学の person-years は壁でなく道。

★同日の指示:

> 形式化作業中に必要なものが出た場合はスケルトン作成とグラフ更新を行ってください。
> 新しい葉を含めたゴール設定を行います

**Why:**

★直前の報告で、[FrdI] の残りを止めている 3 件を「壁」と呼び、
「目標(§2 以降を満点)は現在の在庫では到達できない」と結論していた。
★これは CLAUDE.md の**進め方**そのものに反する——
「工数の大きな塊を壁のように認識しないことを目的としています」と明記されている。

★★**そして実際に、割ったら 2 件の測り違いが出た**(`frdi-three-chains` に詳細):

1. mathlib の探索先パス **`Analysis/Transcendental/` が存在しなかった**
   (実体は `NumberTheory/Transcendental/`)。無いものを探して 0 件を得ていた。
   ★数体上の Siegel の補題は**在庫にあった**。
2. 「この公理からは出ない」の**見る条を間違えていた**。別の条から**証明できた**。

★すなわち「壁」という言葉は、**測定を止める働き**をしていた。

**How to apply:**

- ★到達不能・在庫不足と**報告する前に**、必ず塊を小目標の DAG に割る。
  台帳は JSON(例: `ResearchPaper/frdi-decomposition.json`)、印字は専用の道具
  (例: `tools/frdi-newleaves.mjs`)。**「葉」= 依存が無い、であって易しいという意味ではない**ので
  各節点に大きさを書く。
- ★statement は `Skeleton/` に**型で固定する**(`sorry` は `Skeleton/` では正しい状態)。
  型が付かない分解は、分解できていない。
- ★★`Meta/Claim.lean` の `LeanStatus.absent` に `searched` を書くとき、
  **そのパスが実在するかを確かめる**。存在しないディレクトリを grep した 0 件は測定ではない。
  `ls` で 1 行確かめれば済む。
- ★「割ったら在庫があった/別の条から出た」は繰り返し起きる。
  **1 つの経路で詰まったことを、到達不能と書かない。**

関連: [[frdi-three-chains]] [[leaf-first-with-graph-feedback]] [[leaves-are-measured-not-guessed]]
