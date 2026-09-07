import ABC3.Found.PGC.LubinTateQuotientDescent

/-!
# `𝒪_{K(x)}^H = 𝒪_{K(x₀)}` —— 道 B に残っていた最後の配管

典拠: T. Yoshida, *Local Class Field Theory via Lubin-Tate Theory*
(Ann. Fac. Sci. Toulouse 17-2, 2008; arXiv math/0606108) **Corollary 6.13** と
**Theorem 6.15**(どちらも物理 p.17)。構造化済み原文は
`ResearchPaper/1_Structured/Local Class Field Theory via Lubin-Tate Theory/section-6.html` の
`id="cor-6-13"`(`data-pdf-page="17"`, `data-item="Corollary 6.13"`)と
`id="thm-6-15"`(`data-pdf-page="17"`, `data-item="Theorem 6.15 (Local Kronecker-Weber)"`)。

原文 (Yoshida08 p.17):
> Corollary 6.13. (i) If G ⊲ H, then G^mH/H = (G/H)^m for all m ∈ R[bb]_≥0. (ii) Let K′/K and K′′/K be two Galois extensions with K′K′′/K totally ramified. If Gal(K′/K)^m = Gal(K′′/K)^m = {id} for m ∈ R[bb]_≥0, then Gal(K′K′′/K)^m = {id}. (iii) Let G be abelian. Then |G/G^m| divides (q − 1)q^m−1 for m ∈ Z[bb]_≥0.

原文 (Yoshida08 p.17):
> Theorem 6.15. (Local Kronecker-Weber theorem) Every finite abelian extension of a local field K is a Lubin-Tate extension, i.e. K^LT = K^ab.

## 本ファイルが担当するもの

直前のノード `Found/PGC/LubinTateQuotientDescent.lean` は、B5 の仮定 `hdeg` を完全に外し、
`htriv` を「商の上付き分岐群が `⊥`」(`hbot`)に置き換えたうえで、
**移送の受け口** `upperRamificationGroup_eq_bot_of_equiv` を用意して終わった。
そこに残った穴がただ 1 つ:

> `fixedRing 𝒪_{K(x)} H ≅ 𝒪_{K(x₀)}`(`H = Gal(K(x)/K(x₀))`)の
> **作用を保つ**環同型と、群同型 `G ⧸ H ≃* Gal(K(x₀)/K)`。

★★**本ファイルがそれを埋める。** 中身は「`𝒪_{K(x)}` の `H` 不変元がちょうど `𝒪_{K(x₀)}`」
という 1 行の主張であり、数学ではなく配管である。

## ★設計 —— 「体の話」と「整数環の話」に割った

持ち場の指示どおり、主張を 2 つに割って**抽象核**を 2 本立てた。
★どちらも `PAdicLocalField` も Lubin-Tate も出てこない。

| § | 抽象核 | 語彙 |
|---|---|---|
| §1 | `ringEquivOfRangeEq` —— **共通の環への単射環準同型 2 本の像が一致すれば環同型** | ★**可換環だけ**。分岐・付値・Galois の語が 1 つも出ない |
| §2 | `galQuotientEquiv` —— `Gal(E/k) ⧸ (restrict h).fixingSubgroup ≃* Gal(F/k)` | 体と Galois のみ。★**付値・分岐の語が 1 つも出ない** |

§4 は受け口に流すだけである。立てた主張は 3 本:

| 宣言 | 何を外したか |
|---|---|
| `upperRamificationGroup_quotient_eq_bot_of_adjoin` | ★`hbot` を **B4 の結論そのもの**(`Gal(K(x₀)/K)^m = {id}`)に落とす |
| `le_of_two_quotients_upperRamification_of_adjoin` | ★`E₁` 側の `htriv` が消えた(`E₂` は任意の中間体のまま) |
| `le_of_two_quotients_upperRamification_of_adjoin_pair` | ★★**両側の `htriv` が消えた。素元も `hϖ` も要らない** |

★★**`hbot₀` の型は B4 `upperRamificationGroup_iteratedLubinTatePsi_eq_bot_of_comap` の
結論と字面まで一致している**(`IntermediateField.adjoin K.carrier ({x₀} : Set K.closure)` と
`K.carrier⟮x₀⟯` は同じ項なので、変換無しで `exact` が通ることを別ファイルで確認した)。

★★**2 つの主張(環同型と群同型)が別々の核の系になっている。** 具体層 §3 は
「`𝒪_{K(x)}` の `H` 不変元 ⟺ `K(x₀)` に入る」(`mem_fixedRing_iff_coe_mem`、
Galois 対応 1 回)と「ノルムは中間体をまたいでも変わらない」(`norm_coe_mem`、`rfl`)の
2 つだけを足して、§1 の核に**像の一致**として代入する。

## ★`lean-idioms.md` #69(`adjoinField`/`adjoinIntegers` の境界)をどう避けたか

★**当たらなかった。** `adjoinField` を 1 度も書いていないからである。
本ファイルは終始 `IntermediateField.adjoin K.carrier {x}`(= `K.carrier⟮x⟯`)と
その部分環 `adjoinIntegers K x` だけで書かれている。
★**#69 は「体の構成 `adjoinField` を経由して整数環に降りる」ときに起きる**ので、
最初から整数環側(正確には「`K.closure` への像」側)に寄せれば発生しない。

## ★`lean-idioms.md` #59(2 層をまたぐ `rfl`)をどう避けたか

★**定型 (c)(2 層を `IntermediateField.restrict` + `restrict_algEquiv` に閉じ込める)**を採った。
ただし本ファイルには**もう 1 つの逃げ道**がある:

★★**3 層の部分型をすべて `K.closure` の中の像として潰す。**
`↥(fixedRing ↥(adjoinIntegers K x) H)` は `K.closure` から見て 3 層深いが、
`fixedRingToClosure`(3 つの包含 `RingHom` の合成)で 1 本の環準同型にしてしまえば、
以後は「`K.closure` の元が等しいか」しか問わない。★**`rfl` を層の間で使わない。**
実際に `rfl` で閉じているのは `norm_coe_mem`(中間体 1 層)と
`coe_restrict_algEquiv`(`restrict` の定義が `fieldRange` なので浅い)の 2 本だけで、
どちらも 1 秒以内に返った。

## ★在庫調査の記録(`lean-idioms.md` #158 —— 同じ名前で書いて `already been declared` を出す)

★**11 本の名前を先に `#check` で叩いた。11 本とも `Unknown identifier`** ——
つまり本ファイルの宣言はすべて新規である(直前の実装者は同じ手で 6 本掘り当てたが、
今回は 0 本だった)。叩いた名前:
`norm_coe_mem` / `fixedRingToClosure` / `adjoinIntegersToClosure` /
`mem_fixedRing_iff_coe_mem` / `fixedRingAdjoinEquiv` / `galQuotientEquiv` /
`ringEquivOfRangeEq` / `coe_restrict_algEquiv` / `range_fixedRingToClosure_eq` /
`fixedRingEquivAdjoinIntegers` / `quotientGaloisEquiv`。

★mathlib 側は先に引いた。使ったのは
`IntermediateField.restrict` / `restrict_algEquiv` / `mem_restrict` /
`AlgEquiv.restrictNormalHom` / `restrictNormalHom_surjective` /
`IntermediateField.restrictNormalHom_ker` / `AlgEquiv.restrictNormal_commutes` /
`AlgEquiv.autCongr` / `QuotientGroup.quotientKerEquivOfSurjective` /
`QuotientGroup.quotientMulEquivOfEq` / `IsGalois.fixedField_fixingSubgroup` /
`Normal.of_algEquiv` / `RingHom.rangeRestrict` / `RingEquiv.subringCongr` /
`MulEquiv.irreducible_iff`。★**「mathlib に無い」と書いた箇所は 1 つも無い。**

## ★★#154 を踏まずに `restrictNormalHom` を使えた理由

`lean-idioms.md` #154 は「`AlgEquiv.restrictNormalHom` に `exact` を当てると
200000 heartbeats で落ちる」と言う。★本ファイルは落ちなかった。差は**大きい方の体**である:
#154 が踏んだのは `K₁ := K.closure`(無限次元の代数閉包)だが、
本ファイルは `K₁ := ↥K.carrier⟮x⟯`(**有限次拡大**)で、しかも
`restrictNormalHom` が現れるのは `galQuotientEquiv` の定義 1 か所と
`coe_galQuotientEquiv_apply` の `rfl` 1 か所だけである。
★**`AlgEquiv.restrictNormal_commutes` を先に `have` で取ってから `rw` する**という
#154 の助言どおりの書き方にしてある。

## 逸脱の記録

1. **`L = K`(`n = 1`)に固定**。決定 D29 と B4・B5・直前ノードに合わせた。
   原典の `L = K_n` 一般は扱っていない。
2. **`m` は `k + 1` の形に固定**(`ℕ∞` の切り詰め引き算・除算を書かないため。
   `lean-idioms.md` #102)。§4 の組み立てだけで、§1〜§3 には `m` が
   実数のまま入っている(制限していない)。
3. **`E₁` を `K.carrier⟮x₀⟯` の形に固定**した(§4 の `le_of_two_quotients_...
   _upperRamification_of_adjoin`)。直前ノードの `le_of_two_quotients_upperRamification_of_finrank`
   は `E₁` が任意の中間体でよいが、B4 の結論を流し込むには「`E₁` が単項生成である」ことが要る。
   ★原典の `K^m_x`(Lubin-Tate 拡大)は定義から単項生成なので、これは制限になっていない。
4. **`Gal(K(x₀)/K)` の正規性は仮定である**(`[Normal K.carrier ↥K.carrier⟮x₀⟯]`)。
   原典では `K^m_x/K` が Galois であることに対応する。

## 退化の自己検査

* ★★**環同型が「作用を保つ」ことを落としていない。**
  `map_smul_fixedRingAdjoinEquiv` がそれで、受け口
  `upperRamificationGroup_eq_bot_of_equiv` の `hcompat` にそのまま渡している。
  ★保たない環同型では `i_ϖ` が変わるので移送は**偽**になる。
* ★★**`Algebra 𝒪_K` を保つことも落としていない。**
  `fixedRingAdjoinEquiv_algebraMap` がそれである。
  ★ただし `Algebra 𝒪[K.carrier] ↥(fixedRing …)` は**インスタンスではない**
  (`fixedRingAlgebra` を `letI` で入れる木の流儀。`SMulCommClass` が要る)ので、
  主張は「`c` が `r ∈ 𝒪_K` の像なら `ψ c` も `r` の像」という形にしてある。
  ★これは `letI` 版と同値で、しかも消費側にインスタンスを強要しない。
* ★★**`H` が正規でないと商が作れない。** どこから出しているか:
  §2 の `fixingSubgroup_restrict_normal` が
  「`H = ker(restrictNormalHom)` だから正規」と言っている(核は常に正規)。
  ★その前提は `[Normal k ↥F]`、すなわち**`K(x₀)/K` が正規**であることである。
  ★正規性を落とすと `G ⧸ H` が群にならず、主張が**書けない**。
* ★★**上付きと下付きを取り違えていない。** §4 で移送するのは
  `upperRamificationGroup`(= `G^m`)だけである。下付き `G_n` は 1 度も出てこない。
* ★**同型で移すとき素元を選び直していない。** `ϖ` を勝手な既約元に取ると
  `ψ ϖ = α₀` が言えないので、★**`ϖ := ψ⁻¹(α₀)` と定義**して
  `Irreducible` を `MulEquiv.irreducible_iff` で移している
  (`irreducible_fixedRingAdjoinEquiv_symm`)。
  ★これで「上付き分岐群が素元の取り方に依らない」という**未証明の命題を使わずに済む**。
* ★**`ℕ∞` の切り詰め引き算・除算は 1 つも書いていない。**
-/

namespace ABC3.Found.PGC

open IsLocalRing IsDiscreteValuationRing
open ABC3.Skeleton.PGC
open IntermediateField
open scoped NormedField Valued Classical

/-! ## §1 抽象核 A —— 像が一致する 2 つの単射環準同型は環同型を与える

★**可換環しか出てこない。** 分岐・付値・Galois・局所体の語彙が 1 つも無い。
`lean_check` 0.29 秒・一発で通った。 -/

section RangeEq

variable {A B L : Type*} [CommRing A] [CommRing B] [CommRing L]

/-- ★★**共通の環 `L` へ単射に埋まる 2 つの環が、同じ像を持つなら同型である**。

`A ≃+* f.range = g.range ≃+* B` と 3 段で合成しただけ。
★これが本ファイルの環同型 `fixedRing 𝒪_{K(x)} H ≅ 𝒪_{K(x₀)}` の**唯一の中身**である。
具体層は「像が一致する」ことだけを言えばよくなる。 -/
noncomputable def ringEquivOfRangeEq (f : A →+* L) (g : B →+* L)
    (hf : Function.Injective f) (hg : Function.Injective g)
    (hr : f.range = g.range) : A ≃+* B :=
  ((RingEquiv.ofBijective f.rangeRestrict
      ⟨fun _ _ hab => hf (congrArg Subtype.val hab), f.rangeRestrict_surjective⟩).trans
    (RingEquiv.subringCongr hr)).trans
    (RingEquiv.ofBijective g.rangeRestrict
      ⟨fun _ _ hab => hg (congrArg Subtype.val hab), g.rangeRestrict_surjective⟩).symm

/-- 同型は `L` の中では恒等である —— `g (e a) = f a`。★これが唯一の計算規則。 -/
theorem apply_ringEquivOfRangeEq (f : A →+* L) (g : B →+* L)
    (hf : Function.Injective f) (hg : Function.Injective g)
    (hr : f.range = g.range) (a : A) :
    g (ringEquivOfRangeEq f g hf hg hr a) = f a := by
  set eg := RingEquiv.ofBijective g.rangeRestrict
      ⟨fun _ _ hab => hg (congrArg Subtype.val hab), g.rangeRestrict_surjective⟩ with heg
  have h1 : eg (ringEquivOfRangeEq f g hf hg hr a)
      = RingEquiv.subringCongr hr (f.rangeRestrict a) := by
    rw [ringEquivOfRangeEq]; simp [heg]
  simpa [heg, RingEquiv.ofBijective] using congrArg (Subtype.val : ↥g.range → L) h1

/-- ★★**同型は作用を保つ** —— ただし `A` と `B` の作用は `L` の上の作用ではないので、
「`L` の中で一致する元は、作用させても `L` の中で一致する」(`hact`)を仮定として受け取る。

★`hact` を落とすと結論は**偽**である(作用を保たない環同型が作れてしまう)。 -/
theorem map_smul_ringEquivOfRangeEq {G G' : Type*} [Group G] [Group G']
    [MulSemiringAction G A] [MulSemiringAction G' B]
    (f : A →+* L) (g : B →+* L) (hf : Function.Injective f) (hg : Function.Injective g)
    (hr : f.range = g.range) (e : G ≃* G')
    (hact : ∀ (σ : G) (a : A) (b : B), f a = g b → f (σ • a) = g (e σ • b))
    (σ : G) (a : A) :
    ringEquivOfRangeEq f g hf hg hr (σ • a) = e σ • ringEquivOfRangeEq f g hf hg hr a := by
  apply hg
  rw [apply_ringEquivOfRangeEq]
  exact hact σ a _ (apply_ringEquivOfRangeEq f g hf hg hr a).symm

/-- ★**同型は底の環 `R` からの構造射を保つ** —— `f`・`g` がどちらも `R` の構造射と
可換なら、同型も可換。★`Algebra 𝒪_K` を保つことの抽象版。 -/
theorem map_algebraMap_ringEquivOfRangeEq {R : Type*} [CommRing R] (u : R →+* A) (v : R →+* B)
    (f : A →+* L) (g : B →+* L) (hf : Function.Injective f) (hg : Function.Injective g)
    (hr : f.range = g.range) (huv : ∀ r : R, f (u r) = g (v r)) (r : R) :
    ringEquivOfRangeEq f g hf hg hr (u r) = v r := by
  apply hg
  rw [apply_ringEquivOfRangeEq]
  exact huv r

end RangeEq

/-! ## §2 抽象核 B —— Galois 対応で `Gal(E/k) ⧸ H` を `Gal(F/k)` にする

★**付値・分岐・局所体の語彙が 1 つも出てこない。** 体と Galois だけである。
`lean_check` 0.74 秒(1 度目は `MonoidHom.mk'` の展開だけで失敗、2 度目に通過)。 -/

section GaloisQuotient

variable {k L : Type*} [Field k] [Field L] [Algebra k L] {F E : IntermediateField k L}

/-- `restrict_algEquiv` は `L` の中では恒等である。★`restrict` は包含の `fieldRange` なので
2 層の `rfl` でも浅く、`lean-idioms.md` #59 には触れない。 -/
theorem coe_restrict_algEquiv (h : F ≤ E) (y : ↥F) :
    (((IntermediateField.restrict_algEquiv h y : ↥(IntermediateField.restrict h)) : ↥E) : L)
      = (y : L) := rfl

/-- 逆向きも `L` の中では恒等。 -/
theorem coe_restrict_algEquiv_symm (h : F ≤ E) (v : ↥(IntermediateField.restrict h)) :
    (((IntermediateField.restrict_algEquiv h).symm v : ↥F) : L) = ((v : ↥E) : L) := by
  conv_rhs => rw [← (IntermediateField.restrict_algEquiv h).apply_symm_apply v]
  exact (coe_restrict_algEquiv h _).symm

variable (h : F ≤ E) [Normal k ↥F] [Normal k ↥E]

/-- `F` が正規なら `E/k` の中間体として見た `restrict h` も正規。 -/
instance normal_restrict : Normal k ↥(IntermediateField.restrict h) :=
  Normal.of_algEquiv (IntermediateField.restrict_algEquiv h)

omit [Normal k ↥E] in
/-- `Gal(E/k) → Gal(restrict h/k)` の核は固定部分群(mathlib の
`IntermediateField.restrictNormalHom_ker` そのもの)。 -/
theorem ker_restrictNormalHom_restrict :
    (AlgEquiv.restrictNormalHom (F := k) (K₁ := ↥E) ↥(IntermediateField.restrict h)).ker
      = (IntermediateField.restrict h).fixingSubgroup :=
  IntermediateField.restrictNormalHom_ker (IntermediateField.restrict h)

/-- ★★**`H = Gal(E/F)` は正規部分群である** —— 核は常に正規だから。

★★これが「`H` が正規でないと商が作れない」問題の**出どころ**である。
前提は `[Normal k ↥F]`、すなわち `F/k` が正規であること。 -/
instance fixingSubgroup_restrict_normal :
    ((IntermediateField.restrict h).fixingSubgroup).Normal :=
  ker_restrictNormalHom_restrict h ▸
    (AlgEquiv.restrictNormalHom (F := k) (K₁ := ↥E)
      ↥(IntermediateField.restrict h)).normal_ker

/-- ★★★**Galois 対応の商版** —— `Gal(E/k) ⧸ Gal(E/F) ≃* Gal(F/k)`。

`restrictNormalHom` の全射性(mathlib)＋核(mathlib)＋第一同型定理、
最後に `AlgEquiv.autCongr` で `restrict h` を `F` に戻す。 -/
noncomputable def galQuotientEquiv :
    (↥E ≃ₐ[k] ↥E) ⧸ (IntermediateField.restrict h).fixingSubgroup ≃* (↥F ≃ₐ[k] ↥F) :=
  ((QuotientGroup.quotientMulEquivOfEq (ker_restrictNormalHom_restrict h).symm).trans
    (QuotientGroup.quotientKerEquivOfSurjective _
      (AlgEquiv.restrictNormalHom_surjective (F := k)
        (K₁ := ↥(IntermediateField.restrict h)) ↥E))).trans
    (AlgEquiv.autCongr (IntermediateField.restrict_algEquiv h).symm)

/-- ★★**唯一の計算規則** —— `L` の中では「`σ` を `F` の元に当てる」だけ。

★これが具体層の `hact` を供給する。`AlgEquiv.restrictNormal_commutes` を
先に `have` で取ってから `rw` している(`lean-idioms.md` #154 の助言どおり)。 -/
theorem coe_galQuotientEquiv_apply (σ : ↥E ≃ₐ[k] ↥E) (y : ↥F) :
    ((galQuotientEquiv h (QuotientGroup.mk σ) y : ↥F) : L)
      = ((σ (IntermediateField.inclusion h y) : ↥E) : L) := by
  have hstep : galQuotientEquiv h (QuotientGroup.mk σ) y
      = (IntermediateField.restrict_algEquiv h).symm
        (σ.restrictNormal ↥(IntermediateField.restrict h)
          (IntermediateField.restrict_algEquiv h y)) := rfl
  have hy : (IntermediateField.restrict_algEquiv h y : ↥E) = IntermediateField.inclusion h y :=
    Subtype.ext (coe_restrict_algEquiv h y)
  have hcomm := AlgEquiv.restrictNormal_commutes (F := k) σ ↥(IntermediateField.restrict h)
    (IntermediateField.restrict_algEquiv h y)
  rw [hstep, coe_restrict_algEquiv_symm]
  exact congrArg (fun t : ↥E => (t : L)) (hcomm.trans (congrArg σ hy))

end GaloisQuotient

/-! ## §3 具体層 —— 環同型 `fixedRing 𝒪_{K(x)} H ≃+* 𝒪_{K(x₀)}` -/

section Concrete

variable {p : ℕ} [Fact p.Prime]

/-- ★**中間体をまたいでもノルムは変わらない**(在庫 `norm_mk_of_le` の一般形)。

`↥F` のノルムは `spectralNorm K.carrier K.closure` の制限なので、
`K.closure` の中の元が同じなら同じ値である。★中間体 **1 層**なので `rfl` が速い。 -/
theorem norm_coe_mem (K : PAdicLocalField p) {F F' : IntermediateField K.carrier K.closure}
    (w : ↥F) (hw : (w : K.closure) ∈ F') :
    ‖(⟨(w : K.closure), hw⟩ : ↥F')‖ = ‖w‖ := rfl

/-- `↥(fixedRing 𝒪_{K(x)} H) ↪ K.closure`。★**3 層の部分型を 1 本の環準同型に潰す。**
以後は「`K.closure` の元が等しいか」しか問わないので `lean-idioms.md` #59 に触れない。 -/
noncomputable def fixedRingToClosure (K : PAdicLocalField p) (x : K.closure)
    [FiniteDimensional K.carrier ↥K.carrier⟮x⟯]
    (H : Subgroup (↥K.carrier⟮x⟯ ≃ₐ[K.carrier] ↥K.carrier⟮x⟯)) :
    ↥(fixedRing ↥(adjoinIntegers K x) H) →+* K.closure :=
  (K.carrier⟮x⟯.val.toRingHom).comp
    ((adjoinIntegers K x).subtype.comp (fixedRing ↥(adjoinIntegers K x) H).subtype)

@[simp] theorem fixedRingToClosure_apply (K : PAdicLocalField p) (x : K.closure)
    [FiniteDimensional K.carrier ↥K.carrier⟮x⟯]
    (H : Subgroup (↥K.carrier⟮x⟯ ≃ₐ[K.carrier] ↥K.carrier⟮x⟯))
    (c : ↥(fixedRing ↥(adjoinIntegers K x) H)) :
    fixedRingToClosure K x H c
      = (((c : ↥(adjoinIntegers K x)) : ↥K.carrier⟮x⟯) : K.closure) := rfl

theorem fixedRingToClosure_injective (K : PAdicLocalField p) (x : K.closure)
    [FiniteDimensional K.carrier ↥K.carrier⟮x⟯]
    (H : Subgroup (↥K.carrier⟮x⟯ ≃ₐ[K.carrier] ↥K.carrier⟮x⟯)) :
    Function.Injective (fixedRingToClosure K x H) :=
  fun _ _ hab => Subtype.ext (Subtype.ext (Subtype.ext hab))

/-- `↥𝒪_{K(x₀)} ↪ K.closure`。 -/
noncomputable def adjoinIntegersToClosure (K : PAdicLocalField p) (x₀ : K.closure) :
    ↥(adjoinIntegers K x₀) →+* K.closure :=
  (K.carrier⟮x₀⟯.val.toRingHom).comp (adjoinIntegers K x₀).subtype

@[simp] theorem adjoinIntegersToClosure_apply (K : PAdicLocalField p) (x₀ : K.closure)
    (y : ↥(adjoinIntegers K x₀)) :
    adjoinIntegersToClosure K x₀ y = ((y : ↥K.carrier⟮x₀⟯) : K.closure) := rfl

theorem adjoinIntegersToClosure_injective (K : PAdicLocalField p) (x₀ : K.closure) :
    Function.Injective (adjoinIntegersToClosure K x₀) :=
  fun _ _ hab => Subtype.ext (Subtype.ext hab)

variable (K : PAdicLocalField p) (x x₀ : K.closure)
  [FiniteDimensional K.carrier ↥K.carrier⟮x⟯]

/-- ★★★**核心(体の話)** —— `𝒪_{K(x)}` の元が `H = Gal(K(x)/K(x₀))` で固定されることと、
その元が `K(x₀)` に入ることは同値。

★Galois 対応(`IsGalois.fixedField_fixingSubgroup`)を **1 回**使うだけ。
★「⟸」の向きは Galois を使わない(`mem_fixingSubgroup_iff` だけ)。 -/
theorem mem_fixedRing_iff_coe_mem [IsGalois K.carrier ↥K.carrier⟮x⟯]
    (h : K.carrier⟮x₀⟯ ≤ K.carrier⟮x⟯) (z : ↥(adjoinIntegers K x)) :
    z ∈ fixedRing ↥(adjoinIntegers K x) (IntermediateField.restrict h).fixingSubgroup
      ↔ ((z : ↥K.carrier⟮x⟯) : K.closure) ∈ K.carrier⟮x₀⟯ := by
  have hfix : ∀ ρ : (↥K.carrier⟮x⟯ ≃ₐ[K.carrier] ↥K.carrier⟮x⟯),
      ρ • z = z ↔ ρ (z : ↥K.carrier⟮x⟯) = (z : ↥K.carrier⟮x⟯) := fun ρ => by
    rw [Subtype.ext_iff, coe_smul_adjoinIntegers]
  rw [mem_fixedRing]
  constructor
  · intro hz
    have hmem : (z : ↥K.carrier⟮x⟯)
        ∈ IntermediateField.fixedField (IntermediateField.restrict h).fixingSubgroup := by
      rw [IntermediateField.mem_fixedField_iff]
      intro ρ hρ
      exact (hfix ρ).1 (hz ρ hρ)
    rw [IsGalois.fixedField_fixingSubgroup] at hmem
    exact (IntermediateField.mem_restrict h _).1 hmem
  · intro hz ρ hρ
    exact (hfix ρ).2
      ((IntermediateField.mem_fixingSubgroup_iff _ _).1 hρ _
        ((IntermediateField.mem_restrict h _).2 hz))

/-- ★★★**核心(整数環の話)** —— `K.closure` の中で「`H` 不変な `𝒪_{K(x)}` の元」の集合と
「`𝒪_{K(x₀)}` の元」の集合はまったく同じ。

上の `mem_fixedRing_iff_coe_mem`(体の話)と `norm_coe_mem`(ノルムの話)を
組み合わせただけである。★**これが §1 の抽象核の唯一の入力。** -/
theorem range_fixedRingToClosure_eq [IsGalois K.carrier ↥K.carrier⟮x⟯]
    (h : K.carrier⟮x₀⟯ ≤ K.carrier⟮x⟯) :
    (fixedRingToClosure K x (IntermediateField.restrict h).fixingSubgroup).range
      = (adjoinIntegersToClosure K x₀).range := by
  ext w
  simp only [RingHom.mem_range]
  constructor
  · rintro ⟨c, rfl⟩
    have hmem : (((c : ↥(adjoinIntegers K x)) : ↥K.carrier⟮x⟯) : K.closure) ∈ K.carrier⟮x₀⟯ :=
      (mem_fixedRing_iff_coe_mem K x x₀ h _).1 c.2
    refine ⟨⟨⟨_, hmem⟩, ?_⟩, rfl⟩
    rw [mem_adjoinIntegers_iff', norm_coe_mem K ((c : ↥(adjoinIntegers K x)) : ↥K.carrier⟮x⟯) hmem]
    exact (mem_adjoinIntegers_iff' K x _).1 (c : ↥(adjoinIntegers K x)).2
  · rintro ⟨y, rfl⟩
    have hmemx : (((y : ↥K.carrier⟮x₀⟯)) : K.closure) ∈ K.carrier⟮x⟯ := h (y : ↥K.carrier⟮x₀⟯).2
    have hnorm : (⟨_, hmemx⟩ : ↥K.carrier⟮x⟯) ∈ adjoinIntegers K x := by
      rw [mem_adjoinIntegers_iff', norm_coe_mem K ((y : ↥K.carrier⟮x₀⟯)) hmemx]
      exact (mem_adjoinIntegers_iff' K x₀ _).1 y.2
    refine ⟨⟨⟨_, hnorm⟩, ?_⟩, rfl⟩
    exact (mem_fixedRing_iff_coe_mem K x x₀ h _).2 (y : ↥K.carrier⟮x₀⟯).2

def fixedRingAdjoinEquiv.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 17, item := "Theorem 6.15", sectionId := "thm-6-15" }

/-- ★★★★★★**求める環同型** —— `𝒪_{K(x)}^{Gal(K(x)/K(x₀))} ≅ 𝒪_{K(x₀)}`。

原文 (Yoshida08 p.17):
> Theorem 6.15. (Local Kronecker-Weber theorem) Every finite abelian extension of a local field K is a Lubin-Tate extension, i.e. K^LT = K^ab.

原典 Theorem 6.15 の証明が `Gal(K^m_x/L)` と `Gal(K′K^m_x/L) ⧸ Gal(K′K^m_x/K^m_x)` を
**同一視して**使っている、その同一視である(原典は同一視を明示しない)。

★中身は §1 の `ringEquivOfRangeEq` に `range_fixedRingToClosure_eq` を代入しただけ。 -/
noncomputable def fixedRingAdjoinEquiv [IsGalois K.carrier ↥K.carrier⟮x⟯]
    (h : K.carrier⟮x₀⟯ ≤ K.carrier⟮x⟯) :
    ↥(fixedRing ↥(adjoinIntegers K x) (IntermediateField.restrict h).fixingSubgroup)
      ≃+* ↥(adjoinIntegers K x₀) :=
  ringEquivOfRangeEq _ _
    (fixedRingToClosure_injective K x (IntermediateField.restrict h).fixingSubgroup)
    (adjoinIntegersToClosure_injective K x₀) (range_fixedRingToClosure_eq K x x₀ h)

/-- 同型は `K.closure` の中では恒等。 -/
@[simp] theorem coe_fixedRingAdjoinEquiv [IsGalois K.carrier ↥K.carrier⟮x⟯]
    (h : K.carrier⟮x₀⟯ ≤ K.carrier⟮x⟯)
    (c : ↥(fixedRing ↥(adjoinIntegers K x) (IntermediateField.restrict h).fixingSubgroup)) :
    ((fixedRingAdjoinEquiv K x x₀ h c : ↥K.carrier⟮x₀⟯) : K.closure)
      = (((c : ↥(adjoinIntegers K x)) : ↥K.carrier⟮x⟯) : K.closure) :=
  apply_ringEquivOfRangeEq (fixedRingToClosure K x (IntermediateField.restrict h).fixingSubgroup)
    (adjoinIntegersToClosure K x₀)
    (fixedRingToClosure_injective K x (IntermediateField.restrict h).fixingSubgroup)
    (adjoinIntegersToClosure_injective K x₀) (range_fixedRingToClosure_eq K x x₀ h) c

/-- ★★★**環同型は作用を保つ** —— `G ⧸ H` の作用が `Gal(K(x₀)/K)` の作用に移る。

★★**これを落とすと受け口 `upperRamificationGroup_eq_bot_of_equiv` の `hcompat` が
埋まらず、上付き分岐群の移送が壊れる。**

★`[MulSemiringAction (G ⧸ H) C]` と `hq` は木の流儀(`lean-idioms.md` #165、
在庫 `ramIndex_quotient_mk` と同じ)で**仮定として受け取る**。消費側は
`letI := quotientMulSemiringActionOfTrivial _ smul_fixedRing_eq_self` と
`hq := quotientSMul_mk_fixedRing` を渡せばよい。 -/
theorem map_smul_fixedRingAdjoinEquiv [IsGalois K.carrier ↥K.carrier⟮x⟯]
    [FiniteDimensional K.carrier ↥K.carrier⟮x₀⟯] [Normal K.carrier ↥K.carrier⟮x₀⟯]
    (h : K.carrier⟮x₀⟯ ≤ K.carrier⟮x⟯)
    [MulSemiringAction ((↥K.carrier⟮x⟯ ≃ₐ[K.carrier] ↥K.carrier⟮x⟯)
        ⧸ (IntermediateField.restrict h).fixingSubgroup)
      ↥(fixedRing ↥(adjoinIntegers K x) (IntermediateField.restrict h).fixingSubgroup)]
    (hq : ∀ (σ : ↥K.carrier⟮x⟯ ≃ₐ[K.carrier] ↥K.carrier⟮x⟯)
        (c : ↥(fixedRing ↥(adjoinIntegers K x)
          (IntermediateField.restrict h).fixingSubgroup)),
      (QuotientGroup.mk σ : (↥K.carrier⟮x⟯ ≃ₐ[K.carrier] ↥K.carrier⟮x⟯)
        ⧸ (IntermediateField.restrict h).fixingSubgroup) • c = σ • c)
    (σ : (↥K.carrier⟮x⟯ ≃ₐ[K.carrier] ↥K.carrier⟮x⟯)
      ⧸ (IntermediateField.restrict h).fixingSubgroup)
    (c : ↥(fixedRing ↥(adjoinIntegers K x)
      (IntermediateField.restrict h).fixingSubgroup)) :
    fixedRingAdjoinEquiv K x x₀ h (σ • c)
      = galQuotientEquiv h σ • fixedRingAdjoinEquiv K x x₀ h c := by
  refine map_smul_ringEquivOfRangeEq _ _ _ _ _ (galQuotientEquiv h) ?_ σ c
  intro τ a b hab
  obtain ⟨τ, rfl⟩ := QuotientGroup.mk_surjective τ
  have hb : IntermediateField.inclusion h (b : ↥K.carrier⟮x₀⟯)
      = ((a : ↥(adjoinIntegers K x)) : ↥K.carrier⟮x⟯) :=
    Subtype.ext ((IntermediateField.coe_inclusion h _).trans hab.symm)
  simp only [fixedRingToClosure_apply, adjoinIntegersToClosure_apply, hq, coe_smul_fixedRing,
    coe_smul_adjoinIntegers, coe_galQuotientEquiv_apply, hb]

/-- ★★**環同型は `Algebra 𝒪_K` を保つ** —— `c` が `r ∈ 𝒪_K` の像なら `ψ c` も `r` の像。

★`Algebra 𝒪[K.carrier] ↥(fixedRing …)` は**インスタンスではない**
(`fixedRingAlgebra` を `letI` で入れる木の流儀)ので、
「像である」という形で述べている。★`letI` 版と同値である。
★**これを落とすと分岐指数が移らない。** -/
theorem fixedRingAdjoinEquiv_algebraMap [IsGalois K.carrier ↥K.carrier⟮x⟯]
    (h : K.carrier⟮x₀⟯ ≤ K.carrier⟮x⟯) (r : 𝒪[K.carrier])
    (c : ↥(fixedRing ↥(adjoinIntegers K x) (IntermediateField.restrict h).fixingSubgroup))
    (hc : (c : ↥(adjoinIntegers K x)) = algebraMap 𝒪[K.carrier] ↥(adjoinIntegers K x) r) :
    fixedRingAdjoinEquiv K x x₀ h c = algebraMap 𝒪[K.carrier] ↥(adjoinIntegers K x₀) r := by
  apply adjoinIntegersToClosure_injective K x₀
  rw [adjoinIntegersToClosure_apply, coe_fixedRingAdjoinEquiv, hc]
  rfl

end Concrete

/-! ## §4 受け口へ流す —— `htriv` が消える -/

section Transport

variable {p : ℕ} [Fact p.Prime]

attribute [local instance] isDiscreteValuationRing_adjoinIntegers

variable (K : PAdicLocalField p) (x x₀ : K.closure)
  [FiniteDimensional K.carrier ↥K.carrier⟮x⟯]
  [FiniteDimensional K.carrier ↥K.carrier⟮x₀⟯]

omit [FiniteDimensional K.carrier ↥K.carrier⟮x₀⟯] in
/-- ★**素元は同型の逆像で選ぶ** —— こうすると `ψ ϖ = α₀` が定義から出るので、
「上付き分岐群は素元の取り方に依らない」という未証明の命題を使わずに済む。 -/
theorem irreducible_fixedRingAdjoinEquiv_symm [IsGalois K.carrier ↥K.carrier⟮x⟯]
    (h : K.carrier⟮x₀⟯ ≤ K.carrier⟮x⟯) {α₀ : ↥(adjoinIntegers K x₀)} (hα : Irreducible α₀) :
    Irreducible ((fixedRingAdjoinEquiv K x x₀ h).symm α₀) :=
  (MulEquiv.irreducible_iff (f := (fixedRingAdjoinEquiv K x x₀ h).symm) (x := α₀)).2 hα

def upperRamificationGroup_quotient_eq_bot_of_adjoin.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 17, item := "Corollary 6.13", sectionId := "cor-6-13" }

/-- ★★★★★★★★**`htriv` を消す定理** —— `Gal(K(x₀)/K)^m = {id}` から
`(Gal(K(x)/K) ⧸ Gal(K(x)/K(x₀)))^m = {id}` が出る。

原文 (Yoshida08 p.17):
> Corollary 6.13. (i) If G ⊲ H, then G^mH/H = (G/H)^m for all m ∈ R[bb]_≥0. (ii) Let K′/K and K′′/K be two Galois extensions with K′K′′/K totally ramified. If Gal(K′/K)^m = Gal(K′′/K)^m = {id} for m ∈ R[bb]_≥0, then Gal(K′K′′/K)^m = {id}. (iii) Let G be abelian. Then |G/G^m| divides (q − 1)q^m−1 for m ∈ Z[bb]_≥0.

原典 Corollary 6.13 (ii) の `Gal(K′′/K)^m = {id}` を、木の
「商群 `G ⧸ H` と固定環 `C = 𝒪_{K(x)}^H` の言葉」に翻訳する 1 段である。
★直前ノードが用意した受け口 `upperRamificationGroup_eq_bot_of_equiv` に、
§3 の環同型 `fixedRingAdjoinEquiv` と §2 の群同型 `galQuotientEquiv` を渡すだけ。

★★**素元は `ϖ = ψ⁻¹(α₀)` に固定してある**(退化の自己検査を参照)。 -/
theorem upperRamificationGroup_quotient_eq_bot_of_adjoin
    [IsGalois K.carrier ↥K.carrier⟮x⟯] [Normal K.carrier ↥K.carrier⟮x₀⟯]
    (h : K.carrier⟮x₀⟯ ≤ K.carrier⟮x⟯)
    [Fintype (↥K.carrier⟮x⟯ ≃ₐ[K.carrier] ↥K.carrier⟮x⟯)]
    [Fintype (↥K.carrier⟮x₀⟯ ≃ₐ[K.carrier] ↥K.carrier⟮x₀⟯)]
    [Fintype ((↥K.carrier⟮x⟯ ≃ₐ[K.carrier] ↥K.carrier⟮x⟯)
      ⧸ (IntermediateField.restrict h).fixingSubgroup)]
    [Fintype ↥((IntermediateField.restrict h).fixingSubgroup)]
    [MulSemiringAction ((↥K.carrier⟮x⟯ ≃ₐ[K.carrier] ↥K.carrier⟮x⟯)
        ⧸ (IntermediateField.restrict h).fixingSubgroup)
      ↥(fixedRing ↥(adjoinIntegers K x) (IntermediateField.restrict h).fixingSubgroup)]
    (hq : ∀ (σ : ↥K.carrier⟮x⟯ ≃ₐ[K.carrier] ↥K.carrier⟮x⟯)
        (c : ↥(fixedRing ↥(adjoinIntegers K x)
          (IntermediateField.restrict h).fixingSubgroup)),
      (QuotientGroup.mk σ : (↥K.carrier⟮x⟯ ≃ₐ[K.carrier] ↥K.carrier⟮x⟯)
        ⧸ (IntermediateField.restrict h).fixingSubgroup) • c = σ • c)
    (α₀ : ↥(adjoinIntegers K x₀)) (m : ℝ)
    (hbot : upperRamificationGroup (↥K.carrier⟮x₀⟯ ≃ₐ[K.carrier] ↥K.carrier⟮x₀⟯) α₀ m = ⊥) :
    upperRamificationGroup ((↥K.carrier⟮x⟯ ≃ₐ[K.carrier] ↥K.carrier⟮x⟯)
        ⧸ (IntermediateField.restrict h).fixingSubgroup)
      ((fixedRingAdjoinEquiv K x x₀ h).symm α₀) m = ⊥ := by
  refine upperRamificationGroup_eq_bot_of_equiv (galQuotientEquiv h)
    (fixedRingAdjoinEquiv K x x₀ h)
    (fun σ c => map_smul_fixedRingAdjoinEquiv K x x₀ h hq σ c) _ m ?_
  rw [RingEquiv.apply_symm_apply]
  exact hbot

def le_of_two_quotients_upperRamification_of_adjoin.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 17, item := "Theorem 6.15", sectionId := "thm-6-15" }

/-- ★★★★★★★★★★**Yoshida 2008 Theorem 6.15 の "thus K′ ⊂ Km_x"、`htriv` を完全に外した形**。

原文 (Yoshida08 p.17):
> Theorem 6.15. (Local Kronecker-Weber theorem) Every finite abelian extension of a local field K is a Lubin-Tate extension, i.e. K^LT = K^ab.

直前ノードの `le_of_two_quotients_upperRamification_of_finrank` から、
`E₁` 側の `hbot`(商の上付き分岐群が `⊥`)が

    `hbot₀ : Gal(K(x₀)/K)^m = {id}`

に置き換わっている。★★**これは B4(Proposition 6.14)の結論そのものの形**であり、
`Found/PGC/LubinTateUpperRamificationVanish.lean` の
`upperRamificationGroup_iteratedLubinTatePsi_eq_bot_of_comap` を**そのまま**渡せる。
★素元 `ϖ` も `hϖ` も要らなくなった(`α₀` の既約性から自動で作る)。

★逸脱: `E₁ = K.carrier⟮x₀⟯` の形に固定した(冒頭「逸脱の記録 3」)。
原典の `K^m_x` は単項生成なので制限になっていない。 -/
theorem le_of_two_quotients_upperRamification_of_adjoin
    (ht : IsTotallyRamifiedAdjoin K x)
    [IsGalois K.carrier ↥K.carrier⟮x⟯] [Normal K.carrier ↥K.carrier⟮x₀⟯]
    {α : adjoinIntegers K x}
    (huni : IsLocalRing.maximalIdeal (adjoinIntegers K x) = Ideal.span {α})
    [Fintype (↥K.carrier⟮x⟯ ≃ₐ[K.carrier] ↥K.carrier⟮x⟯)]
    [Fintype (↥K.carrier⟮x₀⟯ ≃ₐ[K.carrier] ↥K.carrier⟮x₀⟯)]
    (habel : ∀ σ τ : (↥K.carrier⟮x⟯ ≃ₐ[K.carrier] ↥K.carrier⟮x⟯), σ * τ = τ * σ)
    {E₂ : IntermediateField K.carrier K.closure}
    (h₁ : K.carrier⟮x₀⟯ ≤ K.carrier⟮x⟯)
    (h₂ : E₂ ≤ K.carrier⟮x⟯)
    (hsup : K.carrier⟮x₀⟯ ⊔ E₂ = K.carrier⟮x⟯)
    [((IntermediateField.restrict h₂).fixingSubgroup).Normal]
    [Fintype ((↥K.carrier⟮x⟯ ≃ₐ[K.carrier] ↥K.carrier⟮x⟯)
      ⧸ (IntermediateField.restrict h₁).fixingSubgroup)]
    [Fintype ((↥K.carrier⟮x⟯ ≃ₐ[K.carrier] ↥K.carrier⟮x⟯)
      ⧸ (IntermediateField.restrict h₂).fixingSubgroup)]
    [MulSemiringAction ((↥K.carrier⟮x⟯ ≃ₐ[K.carrier] ↥K.carrier⟮x⟯)
        ⧸ (IntermediateField.restrict h₁).fixingSubgroup)
      ↥(fixedRing ↥(adjoinIntegers K x) (IntermediateField.restrict h₁).fixingSubgroup)]
    [MulSemiringAction ((↥K.carrier⟮x⟯ ≃ₐ[K.carrier] ↥K.carrier⟮x⟯)
        ⧸ (IntermediateField.restrict h₂).fixingSubgroup)
      ↥(fixedRing ↥(adjoinIntegers K x) (IntermediateField.restrict h₂).fixingSubgroup)]
    (hq : ∀ (σ : ↥K.carrier⟮x⟯ ≃ₐ[K.carrier] ↥K.carrier⟮x⟯)
        (c : ↥(fixedRing ↥(adjoinIntegers K x)
          (IntermediateField.restrict h₁).fixingSubgroup)),
      (QuotientGroup.mk σ : (↥K.carrier⟮x⟯ ≃ₐ[K.carrier] ↥K.carrier⟮x⟯)
        ⧸ (IntermediateField.restrict h₁).fixingSubgroup) • c = σ • c)
    (hq' : ∀ (σ : ↥K.carrier⟮x⟯ ≃ₐ[K.carrier] ↥K.carrier⟮x⟯)
        (c : ↥(fixedRing ↥(adjoinIntegers K x)
          (IntermediateField.restrict h₂).fixingSubgroup)),
      (QuotientGroup.mk σ : (↥K.carrier⟮x⟯ ≃ₐ[K.carrier] ↥K.carrier⟮x⟯)
        ⧸ (IntermediateField.restrict h₂).fixingSubgroup) • c = σ • c)
    {ϖ' : ↥(fixedRing ↥(adjoinIntegers K x) (IntermediateField.restrict h₂).fixingSubgroup)}
    (hϖ' : Irreducible ϖ') (k : ℕ)
    {α₀ : ↥(adjoinIntegers K x₀)}
    (huni₀ : IsLocalRing.maximalIdeal (adjoinIntegers K x₀) = Ideal.span {α₀})
    (hbot₀ : upperRamificationGroup (↥K.carrier⟮x₀⟯ ≃ₐ[K.carrier] ↥K.carrier⟮x₀⟯) α₀
      ((k + 1 : ℕ) : ℝ) = ⊥)
    (hbot' : upperRamificationGroup ((↥K.carrier⟮x⟯ ≃ₐ[K.carrier] ↥K.carrier⟮x⟯)
        ⧸ (IntermediateField.restrict h₂).fixingSubgroup) ϖ' ((k + 1 : ℕ) : ℝ) = ⊥)
    (hdeg : Module.finrank K.carrier ↥K.carrier⟮x₀⟯
      = (Nat.card 𝓀[K.carrier] - 1) * Nat.card 𝓀[K.carrier] ^ k) :
    E₂ ≤ K.carrier⟮x₀⟯ := by
  haveI : Fintype ↥((IntermediateField.restrict h₁).fixingSubgroup) := Fintype.ofFinite _
  have hα₀ : Irreducible α₀ :=
    (IsDiscreteValuationRing.irreducible_iff_uniformizer α₀).2 huni₀
  exact le_of_two_quotients_upperRamification_of_finrank K x ht huni habel h₁ h₂ hsup
    (hq := hq) (hq' := hq')
    (irreducible_fixedRingAdjoinEquiv_symm K x x₀ h₁ hα₀) hϖ' k
    (upperRamificationGroup_quotient_eq_bot_of_adjoin K x x₀ h₁ hq α₀ _ hbot₀) hbot' hdeg

def le_of_two_quotients_upperRamification_of_adjoin_pair.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 17, item := "Theorem 6.15", sectionId := "thm-6-15" }

/-- ★★★★★★★★★★★★**両側とも `htriv` が消えた形** —— B5・直前ノードの
`htriv` / `htriv′`(あるいは `hbot` / `hbot′`)が、どちらも
「その体自身の Galois 群の上付き分岐群が `⊥`」という**原典の字面**になっている。

原文 (Yoshida08 p.17):
> Theorem 6.15. (Local Kronecker-Weber theorem) Every finite abelian extension of a local field K is a Lubin-Tate extension, i.e. K^LT = K^ab.

原典の

    `Gal(K′/K)^m = Gal(K^m_x/K)^m = {id}` ⟹ `K′ ⊂ K^m_x`

そのものである(`L = K`、`m = k+1`)。★★**素元も `hϖ` も 1 つも要らない**
(`huni₀` / `huni₁` から自動で作る)。★★**商群も固定環も表に出てこない。**

★逸脱: `E₁ = K.carrier⟮x₀⟯`、`E₂ = K.carrier⟮x₁⟯` と**どちらも単項生成に固定**した
(冒頭「逸脱の記録 3」)。原典の `K^m_x` は定義から単項生成であり、
`K′` は有限分離拡大なので原始元定理で単項生成になる。★制限になっていない。

★`[((IntermediateField.restrict h₂).fixingSubgroup).Normal]` は**要らなくなった**
——§2 の `fixingSubgroup_restrict_normal` が `[Normal K.carrier ↥K.carrier⟮x₁⟯]` から出す。 -/
theorem le_of_two_quotients_upperRamification_of_adjoin_pair
    (ht : IsTotallyRamifiedAdjoin K x)
    [IsGalois K.carrier ↥K.carrier⟮x⟯] [Normal K.carrier ↥K.carrier⟮x₀⟯]
    {α : adjoinIntegers K x}
    (huni : IsLocalRing.maximalIdeal (adjoinIntegers K x) = Ideal.span {α})
    [Fintype (↥K.carrier⟮x⟯ ≃ₐ[K.carrier] ↥K.carrier⟮x⟯)]
    [Fintype (↥K.carrier⟮x₀⟯ ≃ₐ[K.carrier] ↥K.carrier⟮x₀⟯)]
    (habel : ∀ σ τ : (↥K.carrier⟮x⟯ ≃ₐ[K.carrier] ↥K.carrier⟮x⟯), σ * τ = τ * σ)
    (x₁ : K.closure)
    [FiniteDimensional K.carrier ↥K.carrier⟮x₁⟯] [Normal K.carrier ↥K.carrier⟮x₁⟯]
    [Fintype (↥K.carrier⟮x₁⟯ ≃ₐ[K.carrier] ↥K.carrier⟮x₁⟯)]
    (h₁ : K.carrier⟮x₀⟯ ≤ K.carrier⟮x⟯) (h₂ : K.carrier⟮x₁⟯ ≤ K.carrier⟮x⟯)
    (hsup : K.carrier⟮x₀⟯ ⊔ K.carrier⟮x₁⟯ = K.carrier⟮x⟯)
    [Fintype ((↥K.carrier⟮x⟯ ≃ₐ[K.carrier] ↥K.carrier⟮x⟯)
      ⧸ (IntermediateField.restrict h₁).fixingSubgroup)]
    [Fintype ((↥K.carrier⟮x⟯ ≃ₐ[K.carrier] ↥K.carrier⟮x⟯)
      ⧸ (IntermediateField.restrict h₂).fixingSubgroup)]
    [MulSemiringAction ((↥K.carrier⟮x⟯ ≃ₐ[K.carrier] ↥K.carrier⟮x⟯)
        ⧸ (IntermediateField.restrict h₁).fixingSubgroup)
      ↥(fixedRing ↥(adjoinIntegers K x) (IntermediateField.restrict h₁).fixingSubgroup)]
    [MulSemiringAction ((↥K.carrier⟮x⟯ ≃ₐ[K.carrier] ↥K.carrier⟮x⟯)
        ⧸ (IntermediateField.restrict h₂).fixingSubgroup)
      ↥(fixedRing ↥(adjoinIntegers K x) (IntermediateField.restrict h₂).fixingSubgroup)]
    (hq : ∀ (σ : ↥K.carrier⟮x⟯ ≃ₐ[K.carrier] ↥K.carrier⟮x⟯)
        (c : ↥(fixedRing ↥(adjoinIntegers K x)
          (IntermediateField.restrict h₁).fixingSubgroup)),
      (QuotientGroup.mk σ : (↥K.carrier⟮x⟯ ≃ₐ[K.carrier] ↥K.carrier⟮x⟯)
        ⧸ (IntermediateField.restrict h₁).fixingSubgroup) • c = σ • c)
    (hq' : ∀ (σ : ↥K.carrier⟮x⟯ ≃ₐ[K.carrier] ↥K.carrier⟮x⟯)
        (c : ↥(fixedRing ↥(adjoinIntegers K x)
          (IntermediateField.restrict h₂).fixingSubgroup)),
      (QuotientGroup.mk σ : (↥K.carrier⟮x⟯ ≃ₐ[K.carrier] ↥K.carrier⟮x⟯)
        ⧸ (IntermediateField.restrict h₂).fixingSubgroup) • c = σ • c)
    (k : ℕ)
    {α₀ : ↥(adjoinIntegers K x₀)}
    (huni₀ : IsLocalRing.maximalIdeal (adjoinIntegers K x₀) = Ideal.span {α₀})
    (hbot₀ : upperRamificationGroup (↥K.carrier⟮x₀⟯ ≃ₐ[K.carrier] ↥K.carrier⟮x₀⟯) α₀
      ((k + 1 : ℕ) : ℝ) = ⊥)
    {α₁ : ↥(adjoinIntegers K x₁)}
    (huni₁ : IsLocalRing.maximalIdeal (adjoinIntegers K x₁) = Ideal.span {α₁})
    (hbot₁ : upperRamificationGroup (↥K.carrier⟮x₁⟯ ≃ₐ[K.carrier] ↥K.carrier⟮x₁⟯) α₁
      ((k + 1 : ℕ) : ℝ) = ⊥)
    (hdeg : Module.finrank K.carrier ↥K.carrier⟮x₀⟯
      = (Nat.card 𝓀[K.carrier] - 1) * Nat.card 𝓀[K.carrier] ^ k) :
    K.carrier⟮x₁⟯ ≤ K.carrier⟮x₀⟯ := by
  haveI : Fintype ↥((IntermediateField.restrict h₂).fixingSubgroup) := Fintype.ofFinite _
  have hα₁ : Irreducible α₁ :=
    (IsDiscreteValuationRing.irreducible_iff_uniformizer α₁).2 huni₁
  exact le_of_two_quotients_upperRamification_of_adjoin K x x₀ ht huni habel h₁ h₂ hsup
    (hq := hq) (hq' := hq')
    (irreducible_fixedRingAdjoinEquiv_symm K x x₁ h₂ hα₁) k huni₀ hbot₀
    (upperRamificationGroup_quotient_eq_bot_of_adjoin K x x₁ h₂ hq' α₁ _ hbot₁) hdeg

end Transport

end ABC3.Found.PGC
