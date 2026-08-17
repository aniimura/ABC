---
name: frdi-three-chains
description: [FrdI] の残りを止めている 3 つの塊は「壁」ではなくチェーン(小目標の DAG)。node tools/frdi-newleaves.mjs が葉と層を出す。
metadata:
  type: project
---

★★**2026-08-18 に `frdi-three-walls` を全面的に書き換えたもの。**
旧版は「§2 以降を満点にするのは現在の在庫では到達不能」と結論していた。
**その読み方は撤回した**(CLAUDE.md の姿勢:「工数の山を『壁』と呼ばない。
既知数学の person-years は壁でなく道」)。

## 3 つの塊と、割った先

台帳 `ResearchPaper/frdi-decomposition.json`(22 節点)、印字 `node tools/frdi-newleaves.mjs`。

| チェーン | 奉じる項目 | 下流 | 節点/層/葉 |
|---|---|---|---|
| ★`otricomm` | `Proposition 4.4, (ii)` | **11 件**(最大の律速) | 7 / 5 / 葉 2 は**閉じた** |
| `prol` | `Definition 2.8, (ii)` | 1 件 | 7 / 4 / 葉 4 |
| `sixexp` | `Lemma 6.5, (ii)` | 1 件 | 8 / 6 / 葉 3 |

## ★★割ったら在庫の測定が 2 件覆った(これが「割る」ことの実利)

1. ★**six exponentials**。旧測定「mathlib の `Analysis/Transcendental/` には
   Liouville と `e` しか無い」は、**そのディレクトリが存在しない**(実体は
   `NumberTheory/Transcendental/`)。★実際には
   **`NumberTheory/NumberField/House.lean` に数体上の Siegel の補題**
   (`NumberField.house.exists_ne_zero_int_vec_house_le`)があり、
   これは超越性証明の**算術側の心臓部**である。
2. ★**`Prop 4.4, (ii)` の Ore の四角形**。「(i)(b) `preStepSpan` からは出ない」は
   **見る条の誤り**。効くのは **(iii)(d) `coaPreOverEquiv`**
   (`(𝒞^coa-pre)_A ≃ Order(Φ(A))^opp`、前順序圏 × `Φ(A)` の有向性 ⟹ 余フィルター)。
   ★★**一般の Frobenioid で証明できた**(`Found/FrdI/Prop44Ore.lean`、`sorry` 無し)。

## ★★2026-08-18 の実績 —— 新しい葉 8/8 を同日に閉じた

| チェーン | 閉じた葉 | 実装 |
|---|---|---|
| otricomm | `ore-square` `divgp-hom` `ker-eq-image` | `Found/FrdI/Prop44Ore.lean` `Prop44Ker.lean` |
| prol | `finite-primary` `prol-def` `M-limit` `limit-pi-comm` | `Found/ProL/` 4 ファイル |
| sixexp | `liouville-ineq` `schwarz-many-zeros` `aux-function` `growth` | `Found/SixExp/` 3 ファイル |

★★**次に着手できる葉は 3 個**: `central-ext`(otricomm 層2)/
`M-l-part`(prol 層1)/ `extrapolation`(sixexp 層3)。

## ★`otricomm` の残り(2026-08-18 時点)

`ore-square`(済)・`divgp-hom`(済)・`image-central`(済)・`ker-eq-image`(★済)
→ `central-ext`(中心拡大 ⟹ 交換子が交代双線形形式)
→ ★★**`pairing-vanishes`(その形式が消えること。ここが本体)**
→ `otriBase`(`birat_otriBase_of_comm` で半済)

★★**穴の形が変わった**——「圏論の条件 1 つ」から「**1 つの交代双線形形式が消えるか**」へ。
反例があるならその形式が非零な例として見える。

**How to apply:**

- ★数字を確認するときは `node tools/frdi-newleaves.mjs`(葉と層)と
  `node tools/frdi-blocked.mjs`(下流の件数)を回す。**手で数えない。**
- ★`frdi-blocked.mjs` の出力は**到達不能の証明ではない**。そう読んでいた時期がある。
- ★チェーン待ちがゼロの節は **§1 と §3 だけ**。いちばん早く満点にできるのは **§3**
  (`Proposition 3.2` → `Theorem 3.4`、`Proposition 1.14` 経由)。

関連: [[no-wall-decompose-instead]] [[leaves-are-measured-not-guessed]] [[leaf-first-with-graph-feedback]]
