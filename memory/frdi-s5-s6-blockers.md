---
name: frdi-s5-s6-blockers
description: FrdI §5 と §6 の残りは mathlib の在庫不足 2 群に帰着する（Nikolov–Segal / 因子論・six exponentials）
metadata:
  type: project
---

FrdI の §5（2/7）と §6（0/5）の残りを 2026-08-20 に項目ごと測った。**残りは 2 群の在庫不足に帰着する。**（数え方: `tools/frdi-progress.mjs` は 条なし の `.src`＝命題全体が揃った印だけを数える。個々の段が済んでも命題が揃うまで数字は動かない。）

- **`Proposition 5.6` → Nikolov–Segal**。`M_p ⊆ 𝒪^×(A)` が `Aut_𝒞(A)` の共役で閉じることが要るが、`M_p` は位相（pro-`l` 部分）で定義されるので共役が**連続**でなければ保たれない。位相的有限生成な副有限群で抽象自己同型がすべて連続なのが Nikolov–Segal で、mathlib に `Nikolov` は **0 件**。★これは `Definition 2.8, (i)` の角括弧「[uniquely determined]」（`Def28.lean` が主張しないと明記したもの）を**そのまま消費している**。`Corollary 5.7` はその下流。記録は `Gap/FrdI/Section5.lean` の `Gap_5_6_Mp`。
- **§6 全 5 件 → 因子論と six exponentials**。`WeilDivisor` / `CartierDivisor` / `PrimeDivisor` は mathlib に **0 件**（`Scheme.functionField`・`Morphisms/Proper`・`Geometrically/Integral`・相対正規化 `f.normalization` は**在る**）。`Example 6.1` は `Φ`・`B`・`Div_B` さえ作れば `Theorem 5.2, (ii)`（実装済み）を当てるだけになる。

**それ以外（Prop 5.3 / Cor 5.4 / Prop 5.5）は在庫不足ではなく、残りが 1〜2 の技術的段まで縮んでいる。**

- `Proposition 5.5, (i)` は **2026-08-20 に閉じた**（`Prop55Pf.lean` の `otriPfEquiv`：`𝒪^▷(A)^pf ≃+ 𝒪^▷(A^pf)`。全射性・単射性・積の保存・「自然な関手が定める」まで）。★ただし `A` が Frobenius-trivial かつ Frobenius-normalized の場合＝原文が「immediately」と言う場合。一般の `A` は `Theorem 5.1, (i)` の pre-step の対で降ろす。
- `Corollary 5.4` の `hbirat`（`η ∘ Div^gp = Div^gp ∘ Ψ^birat`）も **2026-08-20 に閉じた**（`Cor54Birat.lean`）。残るのは 1-一意性・rigidity・`Proposition 5.3` の縦の矢印との 1-可換性の 3 条。
- `Proposition 5.3` の「`(𝒞^un-tr)^pf` は model 型」は、**完全化が unit-trivial 性を保つ**こと（`Prop55PfUnit.lean` の `otimes_pfRoot_eq_bot`、Prop 5.5 (i) の同型を単元に制限しただけ）＋ `Theorem 5.1, (iv)` で出る。残りは対象 `(A,n)` への拡張と有理関数単系の同定。
- `Proposition 5.3` は `p55i` 経由の perfection の同値待ち。**行き先**（`Crlf` / `CuntrPf`）と縦の矢印（`untrToSc` / `cToSc`）は作ってある。

★足した仮定は 2 つで、どちらも `index.html` に開示済み：`𝒟` が of FSM-type（`Φ^birat` を monoid on 𝒟 にするため）、`Φ^rlf` の条件 (a)（テンソルの平坦性。perf-factorial なら真と測定済み）。

関連: [[no-wall-decompose-instead]] [[measure-mathlib-before-skeleton]] [[frdi-s1-two-blockers]] [[frdi-three-chains]]

★★2026-08-20 の追記（この日の作業で変わったこと）

- **取りこぼし 3 件の正体は Definition 2.8 / Proposition 1.14 / Proposition 1.6**（`git log --grep=取りこぼし` で特定）。Def 2.8 と Prop 1.14 は閉じた（どちらも逸脱を開示した上での条なし `.src`）。**Prop 1.6 は (ii) を閉じ、残りは 2 点**：(v) metrically trivial の `⇐`（`Gap_1_6_v`）と (vi) `Aut^sub-ample`（sub-automorphism の証人対象を `𝒞′` へ持ち上げる道は `plBk_baseChange` で作れるが、その証人の自己射が同型であることが出ない）。
- **§6 の律速の見積もりを訂正した**。「因子論は mathlib に実質不在で一分野に近い」と書いていたが、**土台は mathlib の部品だけで組めた**。鍵は `IsDiscreteValuationRing.TFAE` の第 4 項「整閉 ∧ 非零素イデアルがちょうど 1 つ」で、これに `isIntegrallyClosed_of_isLocalization` と `IsLocalization.AtPrime.orderIsoOfPrime` を合わせると **任意次元の正規 Noether 整域で DVR** が出る（mathlib は Dedekind の場合しか持たない）。
- その上に `ord`・台の有限性・素因子の型・`div(a)`/`div(f)`・群準同型・有効因子の単系までを `Found/Divisor/HeightOneDVR.lean` に組んだ。**§6 に残るのは Cartier 性・`Q`-Cartier 性、スキーム層への持ち上げ、`V[L]` の正規化と `Φ(L)` の切り出しである。**
- **§6 の底の圏 `𝒟 = B(G)⁰` は閉じた**（`Found/FrdI/Sec6GaloisCat.lean`）。連結・totally epimorphic・of FSM-type。副産物として、`Φ^birat` を monoid on `𝒟` にするために足した (B) の仮定が**幾何の応用では真**である根拠が得られた。

★★★2026-08-20 の追記 2 —— **`Example 6.1` の単系論は全部閉じた**

幾何を落とし、`S` を素因子の型、`Γ ≤ ℤ[S]` を Cartier 因子の群、
`Φ := Γ ∩ ℤ≥0[S]` と置くと、原文の「Observe that …」「one verifies immediately …」
はすべて単系論の定理になる（`Found/Divisor/Cartier*.lean` の 5 ファイル）。

- `isDivisorial_effSub` —— **`Q`-Cartier 性は要らない**。効くのは `Γ` が部分群だけ
- `effSubGpEquiv` —— `Φ^gp ≃ Γ`（ここで初めて `Q`-Cartier が効く）
- `mprec_effSub_iff` —— ★**`a ⪯ b ⇔ supp a ⊆ supp b`**。これが primary の判定を与え、
  `effSubPrimeEquiv : Prime(Φ) ≃ D_L` と台の同定が出る
- `pfEquivNonneg` —— `Φ^pf ≃ ℚ≥0[D_L]`
- `isPerfFactorial_effSub` —— **perf-factorial（11 条全部）**

★鍵は 1 本の定理：**`factorMap ι a p = ι p a`**（`factorMap_iotaAt`）。
原文の `sup(Bound^p_{0}(a))` が `a` の `p` 成分そのものになる。
`≥` 向きは「`a` の `p` 成分だけを取り出した元が実際に `Φ^pf` にあり、
`pCarrierPf p` に入る」ことから出る。

★`realScale`（自分で足した条）はここでは**空虚**である —— `M_p ≃+ ℕ` で、
`ℝ≥0` は可除だが `ℕ` は違うので `ℕ ≃+ ℝ≥0` はない。

★★**`Example 6.1` に残るのは 2 つだけ**：
1. `V[L]`（相対正規化）と `D_L` のスキーム層での構成（`Γ` の実現）
2. `L ↦ Φ(L)` / `L ↦ B(L)` の**関手性**（因子の引き戻し。`monoid on ᴰ` にする）
底の圏は閉じているので、この 2 つが揃えば `Theorem 5.2, (ii)`（実装済み）を当てるだけである。
