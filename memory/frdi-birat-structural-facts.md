---
name: frdi-birat-structural-facts
description: 𝒞^birat を扱うときに効く 3 つの構造的事実。自己射は必ず co-angular、同じ底でも 𝒞^birat で同型とは限らない、群化は行き先が群なら構成不要。
metadata:
  type: project
---

★2026-08-18、`Proposition 4.4` を (i)(ii)(iv) まで閉じる過程で測ったこと。
関連: [[frdi-two-blockers-s3-s4]] [[no-wall-decompose-instead]]

## ★1. 自己射はすべて co-angular

`Definition 1.3, (iii)(b)`(`coAngularOfPreStep`)は
「**co-angular pre-step が 1 本でもあれば、その 2 対象間の射はすべて co-angular**」
と言う。★`𝟙 A` がそれなので、**`A` の自己射はすべて co-angular**。
在庫は `Prop110.lean` の `endo_isCoAngular P F φ`。

★**効き所**: `𝒪^▷(A)` の元は base-identity ＋ linear なので pre-step、
上より co-angular、よって `𝒞^birat` で**同型**になる(`otri_isIso_birat`)。
「写像が立つか」を悩む前にこれを確認すること。

## ★★2. `𝒟` 上の関手性は自分で作らない —— `Proposition 2.2` が供給する

`Φ^birat` を `𝒟` 上の部分関手にするのに「同じ底の対象は `𝒞^birat` で同型」を
証明しようとして遠回りした。★**原典にはその筋は無い**。

FrdI p.83 の `Proposition 4.4, (ii)` 本文が
「the functor `𝒪^×(−)` **on 𝒟** associated to the Frobenioid `𝒞^birat`
[cf. **Proposition 2.2, (ii), (iii)**]」と書いており、
★★**`𝒟` 上の関手性は `Proposition 2.2` が供給する**。
在庫は `Prop22Star.lean` の `otriStar` / `otimes_eq_units` / `dstarMap_otimes_mem`。

★p.85 は「(iii) は (i) の関手の存在と `Proposition 1.5, (ii)` から
immediately に従い、核は `Definition 1.3, (vi)` から出る」と明示している。

★注意した落とし穴: `isIso_of_preStep_of_isGroupLikeObj` は `isotropic` を要求する
(co-angular を `Proposition 1.4, (i)` から取るため)。`𝒞^birat` の対象は
`birat_isIsotropic_iff` より **`𝒞` 側が isotropic なときに限って** isotropic。

★また (iv) の辞書「`𝒞` の射が `𝒞^birat` で同型 ⟺ co-angular pre-step」から
「同じ底でも同型とは限らない」と結論するのは**行き過ぎ**である。
言えるのは「**その射**が同型でない」ことだけで、別の射で同型になりうる。

★副産物: `IsCoAngular` の定義に `γ=𝟙, β=φ, α=𝟙` を入れると
「co-angular な isometric pre-step は同型」が直ちに出る(= `prop_1_4_iii`)。

## ★3. 群化は行き先が群なら構成しない

`𝒪^▷(A)^gp ↪ 𝒪^×(A^birat)` は、行き先が**既に群**なので
**たすき掛けの等式**として書ける: `αβ⁻¹ = γδ⁻¹` ならば `δα = βγ`。

★これで mathlib の `Localization`(可換モノイドの群化)を使わずに済む。
★本プロジェクトの `Gp` は**加法**可換単系用で、`𝒪^▷(A)` は `End A` の
**乗法**部分単系なので、その不一致を回避できるのが大きい。

**Why:** 1 と 3 はどちらも「大きく見えた節点が既存 1 本で片付いた」例で、
2 は逆に「素朴な手順が辞書と矛盾する」例である。
測ってから書くと、3 つとも着手前には見えていなかった。

**How to apply:** `𝒞^birat` で「同型になるか」を問うときは必ず (iv) の辞書に戻る。
群化が出てきたら、まず行き先が既に群でないかを見る。
