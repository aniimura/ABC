import ABC3.Found.PGC.ReciprocityLimitEquivariance
import ABC3.Found.PGC.AbsClosureModules

/-!
# 相互写像の Lubin-Tate データ非依存性(捩れの上)—— pGC Proposition 1.1 が閉じる

典拠: T. Yoshida, *Local Class Field Theory via Lubin-Tate Theory*
(Ann. Fac. Sci. Toulouse 17-2, 2008; arXiv math/0606108) Corollary 4.9(物理 p.9)。
構造化済み原文は `ResearchPaper/1_Structured/Local Class Field Theory via Lubin-Tate Theory/`
の `section-4.html` の `#cor-4-9`。

原文 (Yoshida08 p.9):
> do not depend on f. (We will drop the subscript f and write K^m, ρ_m, K^LT and ρ.)

## ★★★★★何が出たか

直前の節点(`ReciprocityLimitEquivariance.lean`)は `ArtinUnitEquivariance`(壁)を
**名前の付いた仮定ちょうど 1 本** `ReciprocityDatumIndependenceOnTorsion` に還元して
終わっていた。★**本ファイルはその 1 本を無条件で証明する。**

| 宣言 | 内容 |
|---|---|
| `reciprocityDatumIndependenceOnTorsion_holds` | ★★★★★**残っていた穴(無条件)** |
| `artinUnitEquivariance_holds` | ★★★★★**壁が無条件で出た** |
| `cyclotomicCharacter_recoverable_holds` | ★★★★★**pGC Proposition 1.1(無条件)** |

★`#print axioms` はいずれも `[propext, Classical.choice, Quot.sound]` のみである
(公理は増えていない)。

## ★★原典の証明と本ファイルの対応

原典の証明は 4 行である(`.txt` 516–529 行の直読。★**`>` を付けない** ——
`K^m_f` と `Θ^{K̂,×}_{π,π′}` は pdftotext の layout と raw で下付きの位置が食い違うので、
逐語引用には使わない):

    Proof. For f, f′ with linear coefficients π, π′, take θ ∈ Θ and
    [θ] : µ_{f,m} ≅ µ_{f′,m} by Proposition 4.8. Lemma 4.6 shows the two level fields agree.
    If σ(α) = [xπ_j](α) for σ in the Weil group, then
    σ([θ](α)) = [θ](j)[xπ_j](α) = [xπ′_j][θ](α) by Lemma 4.5, hence ρ_{f,m} = ρ_{f′,m}.

| 原典の段 | 本ファイル(または在庫) |
|---|---|
| `θ ∈ Θ^{K̂,×}_{π,π′}` を取る | ★在庫 `exists_arithFrobenius_isCoherent_dworkThetaStep2`(Λ6) |
| `[θ] : µ_{f,m} ≅ µ_{f′,m}` | ★在庫 `exists_bijOn_torsionPoints_of_dworkTheta`(Λ7) |
| `[θ]` が `𝒪_K`-線型(Prop 3.5(iii)) | ★**本ファイル** `subst_comm_of_twisted_intertwine`(抽象核)+ `subst_lubinTateAction_dworkTheta`(具体層) |
| `σ([θ]α) = [θ](σα)`(半線型性、`j = 0` の場合) | ★**本ファイル** `smul_evalAt` / `coe_evalAt_algEquiv` |
| `ρ_{f,m} = ρ_{f′,m}` | ★**本ファイル** `reciprocityMap_eq_of_transport` |
| 極限 `ρ_f = ρ_{f′}` | ★**本ファイル** `reciprocityUnits_eq_of_transport` |

★★**原典の `j ≠ 0`(Weil 群の Frobenius 方向)は本ファイルには出てこない。**
埋めるべき主張が `p` 冪捩れの上に制限されており、捩れ元は自動的に慣性部分
(`j = 0`)に入るからである(§2)。★これが「捩れに制限すると安くなる」ことの中身で、
Lemma 4.5(`uniformizerZ` の 1-コサイクル)を 1 度も消費していない理由でもある。

## ★★★★退化の自己検査 —— 捩れの条件を落とすと偽になる

★★**`ReciprocityDatumIndependenceOnTorsion` から `p^m` 捩れの仮定を落とすと偽である。**
`Gal(K^{ab}/K) ≅ 𝒪_K^× × Ẑ` の分解自体が `π` に依存しており、たとえば
`θ := Art_{π′}(π′)` は `ρ_{f′}(θ) = 1` だが `ρ_f(θ) = π′/π ≠ 1` である。
★捩れ元に限ると `Ẑ` が捩れを持たないことから `θ` は慣性部分に入り(§2)、
慣性の上では素元の選び方に依らない —— これが本ファイルの証明する内容である。

★本ファイルの中核 `reciprocityUnits_eq_of_isUnit_mul` は
「`σ` が `K^{ur}` を各点固定する」という仮定を**明示に**持っており、
これを落とすと(`heqv` の供給が壊れるので)証明が通らない。
★すなわち仮定は飾りではなく、どこで効いているかが 1 箇所に見えている。

## ★★★抽象核 —— `[θ]` の `𝒪_K`-線型性は分岐の言葉を使わない

`subst_comm_of_twisted_intertwine` には**分岐・付値・Galois・Lubin-Tate の語彙が
1 つも出てこない**。一般の可換環 `A` 上の形式冪級数と 1 本の環準同型 `ϕ` だけである:

    `G ∘ θ = θ^ϕ ∘ F`、`F ∘ E = E ∘ F`、`G ∘ E′ = E′ ∘ G`、`E^ϕ = E`、`E′^ϕ = E′`、
    `E` と `E′` の 1 次係数が等しい ⟹ `θ ∘ E = E′ ∘ θ`。

証明は在庫 2 本(`subst_twisted_intertwine_comp` と `powerSeries_uniqueness_twisted`、
どちらも `LubinTateEndoTwisted.lean`)に代入するだけで、**一発で通った**。
★具体層 `subst_lubinTateAction_dworkTheta` は `E := [a]_f`・`E′ := [a]_g`・
`ϕ := σ`(算術 Frobenius)を代入する。`E^ϕ = E` は `map_baseIntHom_fixed`(在庫)、
`F ∘ E = E ∘ F` は `LubinTateAction_functional_equation`(在庫)である。

## ★★段取りとの差分 —— `LubinTateEndoTwisted` は今度こそ要った

★★**本体の見立て(「`ϕ ≠ id` の側なので今度こそ要るかもしれない」)は当たった。**
`powerSeries_uniqueness_twisted`(ねじれ版一意性)と `subst_twisted_intertwine_comp`
(ねじれ絡みの合成)を**どちらも消費している**。★6 波目にして初めて当たった。
★ただし `LubinTateEndoTwisted` そのもの(`[θ]_{f,f′}` の**構成**)は使っていない。
使ったのは**一意性**と**合成**だけである —— `θ` は Λ6(Dwork)が既に作っており、
本ファイルはそれが `𝒪_K` の作用と可換であることだけを一意性から出す。

★**`Θ ≠ ∅`(`surjective_unramGalCompletionUnits_div_self`)は使わなかった。**
Λ6 の `exists_arithFrobenius_isCoherent_dworkThetaStep2` が `θ` を**冪級数として**
直接与えるので、`Θ` の元(係数)を経由する必要が無かった。
★`Θ ≠ ∅` は「係数 `θ` が在る」ことであり、本ファイルが要るのは
「冪級数 `[θ]` が在る」ことである。後者は Λ6 が既に持っていた。

★**素元非依存性の消費のしかた**: 体の一致(`lubinTateClosure_sup_unramifiedClosure_eq_*`)は
**使っていない**。使ったのは Λ7 の**捩れ点の全単射** `exists_bijOn_torsionPoints_of_dworkTheta`
である。★体の一致は「`ρ` の定義域が一致する」ことを言う補題だが、本ファイルの
`reciprocityUnits` は最初から `Gal(K̄/K)` 上の写像なので定義域の問題が起きない。

## ★★新しく要った配管(3 本とも木に無かった)

1. `smul_evalAt` —— `K^{ur}` を各点固定する `σ` は `𝒪_{K̂^{ur}}` 係数の冪級数の
   評価と可換。★`PowerSeries.comp_aeval`(mathlib)に、
   `σ` が `𝒪_{K̂^{ur}}` の像を固定すること(`smul_unramifiedToClosureCompletion`、
   稠密性 + 連続性)から作った `𝒪_{K̂^{ur}}`-代数準同型を渡す。
2. `adjIntToClosureCompletionInt_lubinTateEvalAtTorsionPoint` —— `𝒪_{K(x)}` での
   `aeval` と `𝒪_{ℂ_K}` での `aeval` の一致。★`PowerSeries.hasSum_aeval` を
   連続環準同型で送って `HasSum.unique` で閉じる(`ContinuousSMul 𝒪_K 𝒪_{ℂ_K}` を
   組まずに済む道)。
3. `transport_mem_psi` —— 線型全単射は原始的な捩れ点を原始的な捩れ点へ移す。
   ★`iteratedLubinTateTorsionPoints_sdiff_eq_iteratedLubinTatePsiTorsionPoints`(在庫)と
   `mem_iteratedLubinTateTorsionPoints_iff_forall_mem_maximalIdeal_pow`(在庫、
   ★**`𝔪^m` で書かれているので `π` と `ϖ` の食い違いが起きない**)で出る。

## 逸脱の記録

1. ★**`ϖ = u·π`(`u` は単数)という形に固定した。** 原典は「2 つの素元 `π, π′`」と
   書くが、`Λ6`(Dwork)の出力がこの形でしか無い。一般の `π′` は
   `Ideal.span_singleton_eq_span_singleton` で `u·π` の形に書き直せるので
   (`reciprocityDatumIndependenceOnTorsion_holds` の冒頭 3 行)、一般性は失われていない。
2. ★**`Θ^{K̂}_{π,π′}` の元(係数)を経由していない。** 原典は `θ ∈ Θ` を取って
   `[θ]` を作るが、本ファイルは Λ6 が出す冪級数 `θ` を直接使う。
   両者の関係(`coeff 1 θ ∈ Θ`)は木の `DworkThetaStep2.lean` に既にある。
3. ★**`µ_{f,m}` は `iteratedLubinTateTorsionPoints`(Weierstrass 分解の根集合)で
   表す。** Λ7 と同じ扱いで、両者の一致は
   `aeval_map_baseIntHom_iteratedLubinTate_eq_zero_iff` で確認済み。
4. ★**`W(K^m_f/K)` ではなく `Gal(K̄/K)` の上で述べている。** 埋めるべき
   `ReciprocityDatumIndependenceOnTorsion` がその形だからである。`j` は現れない。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC
open scoped NormedField Valued Classical

variable {p : ℕ} [Fact p.Prime]

def subst_comm_of_twisted_intertwine.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 5, item := "Proposition 3.5", sectionId := "prop-3-5" }

def subst_lubinTateAction_dworkTheta.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 5, item := "Proposition 3.5", sectionId := "prop-3-5" }

def exists_transport_of_dworkTheta.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 9, item := "Corollary 4.9", sectionId := "cor-4-9" }

def reciprocityMap_eq_of_transport.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 9, item := "Corollary 4.9", sectionId := "cor-4-9" }

def reciprocityUnits_eq_of_transport.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 9, item := "Corollary 4.9", sectionId := "cor-4-9" }

def reciprocityUnits_eq_of_isUnit_mul.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 9, item := "Corollary 4.9", sectionId := "cor-4-9" }

def reciprocityDatumIndependenceOnTorsion_holds.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 9, item := "Corollary 4.9", sectionId := "cor-4-9" }

def artinUnitEquivariance_holds.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 3, item := "Proposition 1.1", sectionId := "prop-1-1" }

def cyclotomicCharacter_recoverable_holds.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 3, item := "Proposition 1.1", sectionId := "prop-1-1" }

/-! ## 1. ★★捩れ元は慣性部分に入る

★`Gal(K^{ur}/K) ≅ Ẑ` は捩れを持たないので、`Gal(K^{ab}/K)` の捩れ元は
`K^{ur}` を各点固定する。★これが「捩れに制限すると原典の `j ≠ 0` の段が
まるごと消える」ことの中身である。 -/

/-- `Gal(K^{ur}/K)` は捩れを持たない(`Ẑ` と同型だから)。 -/
theorem unramGal_eq_one_of_pow_eq_one (K : PAdicLocalField p) {n : ℕ} (hn : n ≠ 0)
    {σ : unramGal K} (hσ : σ ^ n = 1) : σ = 1 := by
  have he := (zhatMulEquivUnramGalArith K).symm
  have h1 : (he σ) ^ n = 1 := by rw [← map_pow, hσ, map_one]
  have h2 : he σ = 1 := zhat_eq_one_of_pow_eq_one hn h1
  have := congrArg (zhatMulEquivUnramGalArith K) h2
  simpa using this

/-- ★★**捧れ元の不分岐閉包への制限は自明**。

`Gal(K^{ab}/K)` で `p^m` 捧れなら、`K^{ur} ⊆ K^{ab}` への制限も `p^m` 捧れであり、
`Gal(K^{ur}/K) ≅ Ẑ` は捧れを持たない。 -/
theorem restrictNormalHom_unramifiedClosure_eq_one_of_torsion (K : PAdicLocalField p) (m : ℕ)
    (θ : K.closure ≃ₐ[K.carrier] K.closure)
    (h : (AlgEquiv.restrictNormalHom (abelianClosure K) θ) ^ (p ^ m) = 1) :
    AlgEquiv.restrictNormalHom (unramifiedClosure K) θ = 1 := by
  refine unramGal_eq_one_of_pow_eq_one K (n := p ^ m)
    (pow_ne_zero _ (Fact.out : p.Prime).ne_zero) ?_
  rw [← map_pow]
  have hA : AlgEquiv.restrictNormalHom (abelianClosure K) (θ ^ (p ^ m)) = 1 := by
    rw [map_pow]; exact h
  ext z
  have hyab : (z : K.closure) ∈ abelianClosure K :=
    unramifiedClosure_le_abelianClosure_of_zhat K z.2
  have hfix : (θ ^ (p ^ m)) (z : K.closure) = (z : K.closure) := by
    have hz := AlgEquiv.restrictNormalHom_apply (abelianClosure K) (θ ^ (p ^ m))
      ⟨(z : K.closure), hyab⟩
    rw [hA] at hz
    simpa using hz.symm
  rw [AlgEquiv.restrictNormalHom_apply]
  simpa using hfix

/-- ★捩れ元は `K^{ur}` を各点固定する。 -/
theorem apply_eq_self_of_mem_unramifiedClosure_of_torsion (K : PAdicLocalField p) (m : ℕ)
    (θ : K.closure ≃ₐ[K.carrier] K.closure)
    (h : (AlgEquiv.restrictNormalHom (abelianClosure K) θ) ^ (p ^ m) = 1)
    {y : K.closure} (hy : y ∈ unramifiedClosure K) : θ y = y := by
  have hz := AlgEquiv.restrictNormalHom_apply (unramifiedClosure K) θ ⟨y, hy⟩
  rw [restrictNormalHom_unramifiedClosure_eq_one_of_torsion K m θ h] at hz
  simpa using hz.symm


/-! ## 2. ★★慣性の元は `θ` の評価と可換

`θ` の係数は `𝒪_{K̂^{ur}}` に住むので、`K^{ur}` を各点固定する `σ` は
(連続性と稠密性から)係数を動かさない。したがって `σ(θ(λ)) = θ(σλ)`。
★★原典の「`σ([θ](α)) = [θ]^{(j)}([xπ_j](α))`」の `j = 0` の場合である。 -/

/-- ★`K^{ur}` を各点固定する `σ` は `K̂^{ur}` の像を(`ℂ_K` の中で)各点固定する。

★稠密性(`UniformSpace.Completion.induction_on`)と連続性だけで出る。 -/
theorem smul_unramifiedToClosureCompletion (K : PAdicLocalField p) (σ : K.absGal)
    (hσ : ∀ y ∈ unramifiedClosure K, σ y = y) (w : unramifiedCompletion K) :
    σ • (unramifiedToClosureCompletion K w) = unramifiedToClosureCompletion K w := by
  refine UniformSpace.Completion.induction_on w ?_ ?_
  · exact isClosed_eq
      ((continuous_smul_closureCompletion K σ).comp
        (isometry_unramifiedToClosureCompletion K).continuous)
      (isometry_unramifiedToClosureCompletion K).continuous
  · intro x
    rw [unramifiedToClosureCompletion_coe, smul_coe_closureCompletion]
    exact congrArg _ (hσ (x : K.closure) x.2)

/-- `σ` が `K^{ur}` を各点固定するとき、`ᵒ_{ℂ_K}` への作用は
`ᵒ_{K̂^{ur}}`-代数準同型である。 -/
noncomputable def smulClosureCompletionIntAlgHom (K : PAdicLocalField p) (σ : K.absGal)
    (hσ : ∀ y ∈ unramifiedClosure K, σ y = y) :
    ↥(closureCompletionInt K) →ₐ[↥(unramifiedCompletionInt K)] ↥(closureCompletionInt K) :=
  { MulSemiringAction.toRingHom K.absGal ↥(closureCompletionInt K) σ with
    commutes' := by
      intro r
      apply Subtype.ext
      show ((σ • (algebraMap ↥(unramifiedCompletionInt K) ↥(closureCompletionInt K) r) :
          ↥(closureCompletionInt K)) : closureCompletion K)
        = ((algebraMap ↥(unramifiedCompletionInt K) ↥(closureCompletionInt K) r :
            ↥(closureCompletionInt K)) : closureCompletion K)
      rw [coe_smul_closureCompletionInt, algebraMap_unramifiedInt_closureCompletionInt_coe]
      exact smul_unramifiedToClosureCompletion K σ hσ _ }

/-- 上の代数準同型は作用そのもの。 -/
theorem smulClosureCompletionIntAlgHom_apply (K : PAdicLocalField p) (σ : K.absGal)
    (hσ : ∀ y ∈ unramifiedClosure K, σ y = y) (z : ↥(closureCompletionInt K)) :
    smulClosureCompletionIntAlgHom K σ hσ z = σ • z := rfl

/-- 上の代数準同型は連続。 -/
theorem continuous_smulClosureCompletionIntAlgHom (K : PAdicLocalField p) (σ : K.absGal)
    (hσ : ∀ y ∈ unramifiedClosure K, σ y = y) :
    Continuous (smulClosureCompletionIntAlgHom K σ hσ) := by
  rw [continuous_induced_rng]
  exact (continuous_smul_closureCompletion K σ).comp continuous_subtype_val

/-- ★★`K^{ur}` を各点固定する `σ` は冪級数の評価と可換: `σ(θ(λ)) = θ(σλ)`。 -/
theorem smul_evalAt (K : PAdicLocalField p) (σ : K.absGal)
    (hσ : ∀ y ∈ unramifiedClosure K, σ y = y)
    (θ : PowerSeries ↥(unramifiedCompletionInt K))
    {lam : K.closure} (hlam : ‖lam‖ < 1) (hlam' : ‖σ lam‖ < 1) :
    σ • (evalAt K θ hlam) = evalAt K θ hlam' := by
  haveI := isLinearTopology_closureCompletionInt K
  haveI := continuousSMul_closureCompletionInt K
  have key := AlgHom.congr_fun (PowerSeries.comp_aeval (hasEval_coe_of_norm_lt_one K hlam)
    (continuous_smulClosureCompletionIntAlgHom K σ hσ)) θ
  refine key.trans ?_
  refine aeval_congr_point K _ (hasEval_coe_of_norm_lt_one K hlam') ?_ θ
  apply Subtype.ext
  show ((σ • (⟨(lam : closureCompletion K), _⟩ : ↥(closureCompletionInt K)) :
      ↥(closureCompletionInt K)) : closureCompletion K)
    = (((σ lam : K.closure) : closureCompletion K))
  rw [coe_smul_closureCompletionInt]
  exact smul_coe_closureCompletion K σ lam


/-! ## 3. ★★★抽象核と具体層 —— `[θ]` は `𝒪_K` の作用と可換(Prop 3.5(iii))

原文 (Yoshida08 p.5, Definition 3.3):
> {θ ∈ O_L | θ^ϕ/θ = π′/π}. It is an additive group.

★§3.1 の抽象核には分岐・付値・Galois の語彙が 1 つも出てこない。 -/

/-- ★★★**抽象核** —— ねじれ絡み `θ` は「係数が `ϕ` で動かない自己準同型」と可換。

★分岐・付値・Galois・Lubin-Tate の語彙が 1 つも出てこない。一般の可換環 `A` 上の
形式幂級数と 1 本の環準同型 `ϕ` だけである。

`G ∘ θ = θ^ϕ ∘ F`(`hθ`)、`F ∘ E = E ∘ F`、`G ∘ E′ = E′ ∘ G`、
`E^ϕ = E`、`E′^ϕ = E′`、`E` と `E′` の 1 次係数が等しい —— ならば

  `θ ∘ E = E′ ∘ θ`。

★これが原典 Yoshida Proposition 3.5(iii)(`[θ] ∈ Hom_{ᵒ_L}(F_f, F_{f′})`)の核である。
★退化の自己検査: `hcancel` を落とすと一意性が失われて偽になる
(`powerSeries_uniqueness_twisted` の退化検査を見よ)。 -/
theorem subst_comm_of_twisted_intertwine {A : Type*} [CommRing A] {π π' : A} (ϕ : A →+* A)
    (hcancel : ∀ m : ℕ, 2 ≤ m → ∀ d : A, π' * d = ϕ d * π ^ m → d = 0)
    {F G E E' θ : PowerSeries A}
    (hF0 : PowerSeries.constantCoeff F = 0) (hF1 : PowerSeries.coeff 1 F = π)
    (hG0 : PowerSeries.constantCoeff G = 0) (hG1 : PowerSeries.coeff 1 G = π')
    (hE0 : PowerSeries.constantCoeff E = 0) (hE'0 : PowerSeries.constantCoeff E' = 0)
    (hθ0 : PowerSeries.constantCoeff θ = 0)
    (hEfix : PowerSeries.map ϕ E = E) (hE'fix : PowerSeries.map ϕ E' = E')
    (hEF : PowerSeries.subst E F = PowerSeries.subst F E)
    (hE'G : PowerSeries.subst E' G = PowerSeries.subst G E')
    (hθ : PowerSeries.subst θ G = PowerSeries.subst F (PowerSeries.map ϕ θ))
    (hcoeff : PowerSeries.coeff 1 E = PowerSeries.coeff 1 E') :
    PowerSeries.subst E θ = PowerSeries.subst θ E' := by
  have hEF' : PowerSeries.subst E F = PowerSeries.subst F (PowerSeries.map ϕ E) := by
    rw [hEfix]; exact hEF
  have hE'G' : PowerSeries.subst E' G = PowerSeries.subst G (PowerSeries.map ϕ E') := by
    rw [hE'fix]; exact hE'G
  have hα := subst_twisted_intertwine_comp ϕ hF0 hF0 hG0 hE0 hθ0 hEF' hθ
  have hβ := subst_twisted_intertwine_comp ϕ hF0 hG0 hG0 hθ0 hE'0 hθ hE'G'
  refine powerSeries_uniqueness_twisted ϕ hcancel hF0 hF1
    ((PowerSeries.coeff_zero_eq_constantCoeff_apply G).trans hG0) hG1
    (PowerSeries.constantCoeff_subst_eq_zero hE0 θ hθ0)
    (PowerSeries.constantCoeff_subst_eq_zero hθ0 E' hE'0) ?_ hα hβ
  rw [coeff_one_subst_1var hE0, coeff_one_subst_1var hθ0, hcoeff, mul_comm]

/-- ★★★**具体層** —— Dwork の `θ` は `𝒪_K` の作用と可換: `θ ∘ [a]_f = [a]_g ∘ θ`。 -/
theorem subst_lubinTateAction_dworkTheta (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal 𝒪[K.carrier]) 𝒪[K.carrier]]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField 𝒪[K.carrier]) pp]
    [Fintype (IsLocalRing.ResidueField 𝒪[K.carrier])]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField 𝒪[K.carrier]) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal 𝒪[K.carrier] = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries 𝒪[K.carrier]) (hf0 : PowerSeries.coeff 0 f = 0)
    (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue 𝒪[K.carrier]) f = PowerSeries.X ^ (pp ^ ff))
    {ϖ : 𝒪[K.carrier]} (hϖmax : IsLocalRing.maximalIdeal 𝒪[K.carrier] = Ideal.span {ϖ})
    (hϖne0 : ϖ ≠ 0)
    (g : PowerSeries 𝒪[K.carrier]) (hg0 : PowerSeries.coeff 0 g = 0)
    (hg1 : PowerSeries.coeff 1 g = ϖ)
    (hg : PowerSeries.map (IsLocalRing.residue 𝒪[K.carrier]) g = PowerSeries.X ^ (pp ^ ff))
    (σ : unramGal K) (θ : PowerSeries ↥(unramifiedCompletionInt K))
    (hθ0 : PowerSeries.constantCoeff θ = 0)
    (hint : PowerSeries.subst (PowerSeries.map (baseIntHom K) f)
        (PowerSeries.map (unramGalCompletionIntHom K σ) θ)
      = PowerSeries.subst θ (PowerSeries.map (baseIntHom K) g))
    (a : 𝒪[K.carrier]) :
    PowerSeries.subst
        (PowerSeries.map (baseIntHom K) (LubinTateAction hq hπmax f hf0 hf1 hf a)) θ
      = PowerSeries.subst θ
        (PowerSeries.map (baseIntHom K) (LubinTateAction hq hϖmax g hg0 hg1 hg a)) := by
  haveI := isPrincipalIdealRing_unramifiedCompletionInt K
  have hπ0 : (π : K.carrier) ≠ 0 := by exact_mod_cast hπne0
  have hϖ0 : (ϖ : K.carrier) ≠ 0 := by exact_mod_cast hϖne0
  have hmaxπ : IsLocalRing.maximalIdeal ↥(unramifiedCompletionInt K)
      = Ideal.span {baseIntHom K π} := by
    rw [baseIntHom_eq_uniformizerCompletionInt]
    exact maximalIdeal_unramifiedCompletionInt_eq_span K hπ0 hπmax
  have hmaxϖ : IsLocalRing.maximalIdeal ↥(unramifiedCompletionInt K)
      = Ideal.span {baseIntHom K ϖ} := by
    rw [baseIntHom_eq_uniformizerCompletionInt]
    exact maximalIdeal_unramifiedCompletionInt_eq_span K hϖ0 hϖmax
  have hcancel := semilinear_cancel_of_isNoetherian hmaxϖ
    (by rw [baseIntHom_eq_uniformizerCompletionInt]
        exact uniformizerCompletionInt_ne_zero K hϖ0)
    (by rw [hmaxπ]; exact Ideal.mem_span_singleton_self _)
    (unramGalCompletionIntHom K σ)
    (maximalIdeal_map_mem_of_ringEquiv (unramGalCompletionInt K σ))
  have hfe := congrArg (PowerSeries.map (baseIntHom K))
    (LubinTateAction_functional_equation hq hπmax f hf0 hf1 hf a)
  have hge := congrArg (PowerSeries.map (baseIntHom K))
    (LubinTateAction_functional_equation hq hϖmax g hg0 hg1 hg a)
  rw [map_subst_powerSeries _
      (PowerSeries.HasSubst.of_constantCoeff_zero'
        (constantCoeff_LubinTateAction hq hπmax f hf0 hf1 hf a)),
    map_subst_powerSeries _
      (PowerSeries.HasSubst.of_constantCoeff_zero'
        ((PowerSeries.coeff_zero_eq_constantCoeff_apply f).symm.trans hf0))] at hfe
  rw [map_subst_powerSeries _
      (PowerSeries.HasSubst.of_constantCoeff_zero'
        (constantCoeff_LubinTateAction hq hϖmax g hg0 hg1 hg a)),
    map_subst_powerSeries _
      (PowerSeries.HasSubst.of_constantCoeff_zero'
        ((PowerSeries.coeff_zero_eq_constantCoeff_apply g).symm.trans hg0))] at hge
  refine subst_comm_of_twisted_intertwine (unramGalCompletionIntHom K σ) hcancel
    (constantCoeff_map_eq_zero _ ((PowerSeries.coeff_zero_eq_constantCoeff_apply f).symm.trans hf0))
    (by rw [PowerSeries.coeff_map, hf1])
    (constantCoeff_map_eq_zero _ ((PowerSeries.coeff_zero_eq_constantCoeff_apply g).symm.trans hg0))
    (by rw [PowerSeries.coeff_map, hg1])
    (constantCoeff_map_eq_zero _ (constantCoeff_LubinTateAction hq hπmax f hf0 hf1 hf a))
    (constantCoeff_map_eq_zero _ (constantCoeff_LubinTateAction hq hϖmax g hg0 hg1 hg a))
    hθ0
    (map_baseIntHom_fixed K σ _) (map_baseIntHom_fixed K σ _)
    hfe hge hint.symm ?_
  rw [PowerSeries.coeff_map, PowerSeries.coeff_map,
    coeff_one_LubinTateAction hq hπmax f hf0 hf1 hf a,
    coeff_one_LubinTateAction hq hϖmax g hg0 hg1 hg a]


/-! ## 4. ★★転送 —— `𝒪_{K(x)}` での評価と `𝒪_{ℂ_K}` での評価をつなぐ -/

/-- `ᵒ_{K(x)} → ᵒ_{K^{al}}`。 -/
noncomputable def adjIntToAbsClosureInt (K : PAdicLocalField p) (x : K.closure) :
    adjoinIntegers K x →+* ↥(absClosureInt K) where
  toFun z := ⟨((z : ↥(IntermediateField.adjoin K.carrier ({x} : Set K.closure))) : K.closure),
    (mem_absClosureInt K _).mpr z.2⟩
  map_one' := rfl
  map_mul' := fun _ _ => rfl
  map_zero' := rfl
  map_add' := fun _ _ => rfl

/-- `ᵒ_{K(x)} → ᵒ_{ℂ_K}`(`ᵒ_K`-代数準同型)。 -/
noncomputable def adjIntToClosureCompletionInt (K : PAdicLocalField p) (x : K.closure) :
    adjoinIntegers K x →ₐ[𝒪[K.carrier]] ↥(closureCompletionInt K) :=
  { (absClosureIntToClosureCompletionInt K).comp (adjIntToAbsClosureInt K x) with
    commutes' := by
      intro r
      apply Subtype.ext
      show ((absClosureIntToClosureCompletionInt K) (adjIntToAbsClosureInt K x
          (algebraMap ↥𝒪[K.carrier] ↥(adjoinIntegers K x) r)) : closureCompletion K)
        = ((algebraMap ↥𝒪[K.carrier] ↥(closureCompletionInt K) r : ↥(closureCompletionInt K))
            : closureCompletion K)
      rw [coe_absClosureIntToClosureCompletionInt, coe_algebraMap_base_closureCompletionInt]
      rfl }

/-- 上の写像の値は `ι_{al}` の像そのもの。 -/
@[simp] theorem coe_adjIntToClosureCompletionInt (K : PAdicLocalField p) (x : K.closure)
    (z : adjoinIntegers K x) :
    ((adjIntToClosureCompletionInt K x z : ↥(closureCompletionInt K)) : closureCompletion K)
      = closureCompletionCoe K
        ((z : ↥(IntermediateField.adjoin K.carrier ({x} : Set K.closure))) : K.closure) := rfl

/-- 上の写像は連続(等長の合成)。 -/
theorem continuous_adjIntToClosureCompletionInt (K : PAdicLocalField p) (x : K.closure) :
    Continuous (adjIntToClosureCompletionInt K x) := by
  rw [continuous_induced_rng]
  exact (isometry_closureCompletionCoe K).continuous.comp
    (Isometry.continuous fun _ => congrFun rfl)

/-- スカラー塔: `(bh a) • w = a • w`。 -/
theorem baseIntHom_smul_closureCompletionInt (K : PAdicLocalField p) (a : 𝒪[K.carrier])
    (w : ↥(closureCompletionInt K)) : (baseIntHom K a) • w = a • w := by
  rw [← IsScalarTower.algebraMap_smul ↥(unramifiedCompletionInt K) a w]; rfl

/-- ★★**転送** —— `ᵒ_K` 係数の幂級数を `Λ_n` の点で評価した値は、
`ᵒ_{ℂ_K}` の中の評価と一致する。

★`PowerSeries.hasSum_aeval` を連続環準同型で送って `HasSum.unique` で閉じる。
★**`ContinuousSMul ᵒ_K ᵒ_{ℂ_K}` を組まずに済む**のがこの道の利点。 -/
theorem adjIntToClosureCompletionInt_lubinTateEvalAtTorsionPoint (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal 𝒪[K.carrier]) 𝒪[K.carrier]]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField 𝒪[K.carrier]) pp]
    [Fintype (IsLocalRing.ResidueField 𝒪[K.carrier])]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField 𝒪[K.carrier]) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal 𝒪[K.carrier] = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries 𝒪[K.carrier]) (hf0 : PowerSeries.coeff 0 f = 0)
    (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue 𝒪[K.carrier]) f = PowerSeries.X ^ (pp ^ ff))
    (n : ℕ) (x : K.closure)
    (hx : x ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n)
    (hmem : x ∈ IntermediateField.adjoin K.carrier ({x} : Set K.closure))
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    (hnorm : ‖x‖ < 1) (F : PowerSeries 𝒪[K.carrier]) :
    adjIntToClosureCompletionInt K x
        (lubinTateEvalAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf n x hx hmem F)
      = evalAt K (PowerSeries.map (baseIntHom K) F) hnorm := by
  haveI := completeSpace_adjoinIntegers K x
  haveI := isLinearTopology_adjoinIntegers K x
  haveI := continuousSMul_adjoinIntegers K x
  haveI := isLinearTopology_closureCompletionInt K
  haveI := continuousSMul_closureCompletionInt K
  have hs1 := PowerSeries.hasSum_aeval
    (hasEval_mem_adjoinIntegers_of_mem_iteratedLubinTateTorsionPoints
      K hq hπmax hπne0 f hf0 hf1 hf n x hx hmem) F
  have hs2 := hs1.map (adjIntToClosureCompletionInt K x).toRingHom.toAddMonoidHom
    (continuous_adjIntToClosureCompletionInt K x)
  have hs3 := PowerSeries.hasSum_aeval (hasEval_coe_of_norm_lt_one K hnorm)
    (PowerSeries.map (baseIntHom K) F)
  have hgoal : evalAt K (PowerSeries.map (baseIntHom K) F) hnorm
      = PowerSeries.aeval (hasEval_coe_of_norm_lt_one K hnorm)
          (PowerSeries.map (baseIntHom K) F) := rfl
  rw [hgoal]
  refine hs2.unique ?_
  convert hs3 using 2 with d
  · rfl
  · rw [PowerSeries.coeff_map, baseIntHom_smul_closureCompletionInt]
    show (adjIntToClosureCompletionInt K x) (PowerSeries.coeff d F • _ ^ d) = _
    rw [map_smul, map_pow]
    rfl


/-! ## 5. `K^{al}` の中で見た Lubin-Tate 作用 -/


/-- `K^{al}` の中で見た `[a]_f · λ`。 -/
noncomputable def ltActPt (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal 𝒪[K.carrier]) 𝒪[K.carrier]]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField 𝒪[K.carrier]) pp]
    [Fintype (IsLocalRing.ResidueField 𝒪[K.carrier])]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField 𝒪[K.carrier]) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal 𝒪[K.carrier] = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries 𝒪[K.carrier]) (hf0 : PowerSeries.coeff 0 f = 0)
    (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue 𝒪[K.carrier]) f = PowerSeries.X ^ (pp ^ ff))
    (n : ℕ) (lam : K.closure)
    (hlam : lam ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n)
    (a : 𝒪[K.carrier]) : K.closure :=
  ((lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf n lam hlam
      (IntermediateField.mem_adjoin_simple_self K.carrier lam) a :
    ↥(IntermediateField.adjoin K.carrier ({lam} : Set K.closure))) : K.closure)

/-- `[a]_f · λ` もまた `Λ_n` の元。 -/
theorem ltActPt_mem (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal 𝒪[K.carrier]) 𝒪[K.carrier]]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField 𝒪[K.carrier]) pp]
    [Fintype (IsLocalRing.ResidueField 𝒪[K.carrier])]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField 𝒪[K.carrier]) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal 𝒪[K.carrier] = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries 𝒪[K.carrier]) (hf0 : PowerSeries.coeff 0 f = 0)
    (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue 𝒪[K.carrier]) f = PowerSeries.X ^ (pp ^ ff))
    (n : ℕ) (lam : K.closure)
    (hlam : lam ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n)
    (a : 𝒪[K.carrier]) :
    ltActPt K hq hπmax hπne0 f hf0 hf1 hf n lam hlam a
      ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n :=
  lubinTateActionAtTorsionPoint_mem K hq hπmax hπne0 f hf0 hf1 hf n lam hlam _ a

/-- `[a]_f · λ` の `ℂ_K` での姿は `ᵒ_{ℂ_K}` での評価である。 -/
theorem coe_evalAt_map_lubinTateAction (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal 𝒪[K.carrier]) 𝒪[K.carrier]]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField 𝒪[K.carrier]) pp]
    [Fintype (IsLocalRing.ResidueField 𝒪[K.carrier])]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField 𝒪[K.carrier]) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal 𝒪[K.carrier] = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries 𝒪[K.carrier]) (hf0 : PowerSeries.coeff 0 f = 0)
    (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue 𝒪[K.carrier]) f = PowerSeries.X ^ (pp ^ ff))
    (n : ℕ) (lam : K.closure)
    (hlam : lam ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n)
    (hnorm : ‖lam‖ < 1) (a : 𝒪[K.carrier]) :
    ((evalAt K (PowerSeries.map (baseIntHom K) (LubinTateAction hq hπmax f hf0 hf1 hf a)) hnorm :
        ↥(closureCompletionInt K)) : closureCompletion K)
      = closureCompletionCoe K (ltActPt K hq hπmax hπne0 f hf0 hf1 hf n lam hlam a) := by
  rw [← adjIntToClosureCompletionInt_lubinTateEvalAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf
    n lam hlam (IntermediateField.mem_adjoin_simple_self K.carrier lam) hnorm _]
  rfl


/-! ## 6. ★★★捩れ点の上での `𝒪_K`-線型性 -/

/-- ★★★**`θ([a]_f λ) = [a]_g(θ λ)`** —— 捩れ点の上での `𝒪_K`-線型性。 -/
theorem coe_evalAt_ltActPt (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal 𝒪[K.carrier]) 𝒪[K.carrier]]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField 𝒪[K.carrier]) pp]
    [Fintype (IsLocalRing.ResidueField 𝒪[K.carrier])]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField 𝒪[K.carrier]) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal 𝒪[K.carrier] = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries 𝒪[K.carrier]) (hf0 : PowerSeries.coeff 0 f = 0)
    (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue 𝒪[K.carrier]) f = PowerSeries.X ^ (pp ^ ff))
    {ϖ : 𝒪[K.carrier]} (hϖmax : IsLocalRing.maximalIdeal 𝒪[K.carrier] = Ideal.span {ϖ})
    (hϖne0 : ϖ ≠ 0)
    (g : PowerSeries 𝒪[K.carrier]) (hg0 : PowerSeries.coeff 0 g = 0)
    (hg1 : PowerSeries.coeff 1 g = ϖ)
    (hg : PowerSeries.map (IsLocalRing.residue 𝒪[K.carrier]) g = PowerSeries.X ^ (pp ^ ff))
    (σ : unramGal K) (θ : PowerSeries ↥(unramifiedCompletionInt K))
    (hθ0 : PowerSeries.constantCoeff θ = 0)
    (hint : PowerSeries.subst (PowerSeries.map (baseIntHom K) f)
        (PowerSeries.map (unramGalCompletionIntHom K σ) θ)
      = PowerSeries.subst θ (PowerSeries.map (baseIntHom K) g))
    (n : ℕ) (lam : K.closure)
    (hlam : lam ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n)
    (hlamnorm : ‖lam‖ < 1) (a : 𝒪[K.carrier])
    (hμnorm : ‖ltActPt K hq hπmax hπne0 f hf0 hf1 hf n lam hlam a‖ < 1)
    (elam : K.closure)
    (helam : elam ∈ iteratedLubinTateTorsionPoints K hq hϖmax hϖne0 g hg0 hg1 hg n)
    (helamnorm : ‖elam‖ < 1)
    (hev : ((evalAt K θ hlamnorm : ↥(closureCompletionInt K)) : closureCompletion K)
      = closureCompletionCoe K elam) :
    ((evalAt K θ hμnorm : ↥(closureCompletionInt K)) : closureCompletion K)
      = closureCompletionCoe K (ltActPt K hq hϖmax hϖne0 g hg0 hg1 hg n elam helam a) := by
  haveI := isLinearTopology_closureCompletionInt K
  haveI := continuousSMul_closureCompletionInt K
  have hEa0 : PowerSeries.constantCoeff
      (PowerSeries.map (baseIntHom K) (LubinTateAction hq hπmax f hf0 hf1 hf a)) = 0 :=
    constantCoeff_map_eq_zero _ (constantCoeff_LubinTateAction hq hπmax f hf0 hf1 hf a)
  have hEa'0 : PowerSeries.constantCoeff
      (PowerSeries.map (baseIntHom K) (LubinTateAction hq hϖmax g hg0 hg1 hg a)) = 0 :=
    constantCoeff_map_eq_zero _ (constantCoeff_LubinTateAction hq hϖmax g hg0 hg1 hg a)
  have hlin := subst_lubinTateAction_dworkTheta K hq hπmax hπne0 f hf0 hf1 hf
    hϖmax hϖne0 g hg0 hg1 hg σ θ hθ0 hint a
  have hz := hasEval_coe_of_norm_lt_one K hlamnorm
  have hv := hasEval_evalAt K hlamnorm _ hEa0
  have hw := hasEval_evalAt K hlamnorm θ hθ0
  have c1 : PowerSeries.aeval hz (PowerSeries.subst
        (PowerSeries.map (baseIntHom K) (LubinTateAction hq hπmax f hf0 hf1 hf a)) θ)
      = PowerSeries.aeval hv θ :=
    aeval_subst_eq_aeval_aeval (PowerSeries.HasSubst.of_constantCoeff_zero' hEa0) hEa0 hz hv rfl
  have c2 : PowerSeries.aeval hz (PowerSeries.subst θ
        (PowerSeries.map (baseIntHom K) (LubinTateAction hq hϖmax g hg0 hg1 hg a)))
      = PowerSeries.aeval hw _ :=
    aeval_subst_eq_aeval_aeval (PowerSeries.HasSubst.of_constantCoeff_zero' hθ0) hθ0 hz hw rfl
  have key : PowerSeries.aeval hv θ
      = PowerSeries.aeval hw
        (PowerSeries.map (baseIntHom K) (LubinTateAction hq hϖmax g hg0 hg1 hg a)) := by
    rw [← c1, hlin, c2]
  have h1 : PowerSeries.aeval hv θ = evalAt K θ hμnorm := by
    refine aeval_congr_point K hv (hasEval_coe_of_norm_lt_one K hμnorm) ?_ θ
    apply Subtype.ext
    exact coe_evalAt_map_lubinTateAction K hq hπmax hπne0 f hf0 hf1 hf n lam hlam hlamnorm a
  have h2 : PowerSeries.aeval hw
        (PowerSeries.map (baseIntHom K) (LubinTateAction hq hϖmax g hg0 hg1 hg a))
      = evalAt K (PowerSeries.map (baseIntHom K) (LubinTateAction hq hϖmax g hg0 hg1 hg a))
          helamnorm := by
    refine aeval_congr_point K hw (hasEval_coe_of_norm_lt_one K helamnorm) ?_ _
    apply Subtype.ext
    exact hev
  rw [← h1, key, h2]
  exact coe_evalAt_map_lubinTateAction K hq hϖmax hϖne0 g hg0 hg1 hg n elam helam helamnorm a


/-! ## 7. ★★★★★輸送の完成形

原文 (Yoshida08 p.9):
> do not depend on f. (We will drop the subscript f and write K^m, ρ_m, K^LT and ρ.)

Λ7 の全単射 `Λ_{f,n} ≃ Λ_{g,n}` が、慣性の Galois 作用と可換で、しかも
`𝒪_K`-線型であることを 1 つにまとめる。 -/

/-- ★★`K^{ur}` を固定する `τ` は `θ` の値と可換(`K^{al}` の元としての形)。 -/
theorem coe_evalAt_algEquiv (K : PAdicLocalField p) (τ : K.absGal)
    (hτ : ∀ y ∈ unramifiedClosure K, τ y = y)
    (θ : PowerSeries ↥(unramifiedCompletionInt K))
    {lam : K.closure} (hlamnorm : ‖lam‖ < 1) (hτnorm : ‖τ lam‖ < 1)
    (elam : K.closure)
    (hev : ((evalAt K θ hlamnorm : ↥(closureCompletionInt K)) : closureCompletion K)
      = closureCompletionCoe K elam) :
    ((evalAt K θ hτnorm : ↥(closureCompletionInt K)) : closureCompletion K)
      = closureCompletionCoe K (τ elam) := by
  rw [← smul_evalAt K τ hτ θ hlamnorm hτnorm, coe_smul_closureCompletionInt, hev]
  exact smul_coe_closureCompletion K τ elam

/-- ★★★★★**輸送の完成形** —— Dwork の `θ` が誘導する `Λ_{f,n} ≃ Λ_{g,n}` は
不分岐部分の Galois 作用と可換で、しかも `𝒪_K`-線型である。 -/
theorem exists_transport_of_dworkTheta (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal 𝒪[K.carrier]) 𝒪[K.carrier]]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField 𝒪[K.carrier]) pp]
    [Fintype (IsLocalRing.ResidueField 𝒪[K.carrier])]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField 𝒪[K.carrier]) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal 𝒪[K.carrier] = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries 𝒪[K.carrier]) (hf0 : PowerSeries.coeff 0 f = 0)
    (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue 𝒪[K.carrier]) f = PowerSeries.X ^ (pp ^ ff))
    {ϖ : 𝒪[K.carrier]} (hϖmax : IsLocalRing.maximalIdeal 𝒪[K.carrier] = Ideal.span {ϖ})
    (hϖne0 : ϖ ≠ 0)
    (g : PowerSeries 𝒪[K.carrier]) (hg0 : PowerSeries.coeff 0 g = 0)
    (hg1 : PowerSeries.coeff 1 g = ϖ)
    (hg : PowerSeries.map (IsLocalRing.residue 𝒪[K.carrier]) g = PowerSeries.X ^ (pp ^ ff))
    (σ : unramGal K) (θ θ' : PowerSeries ↥(unramifiedCompletionInt K))
    (hθ0 : PowerSeries.constantCoeff θ = 0)
    (hθ'θ : PowerSeries.subst θ θ' = PowerSeries.X)
    (hint : PowerSeries.subst (PowerSeries.map (baseIntHom K) f)
        (PowerSeries.map (unramGalCompletionIntHom K σ) θ)
      = PowerSeries.subst θ (PowerSeries.map (baseIntHom K) g))
    (n : ℕ) :
    ∃ e : K.closure → K.closure,
      Set.BijOn e ↑(iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n)
        ↑(iteratedLubinTateTorsionPoints K hq hϖmax hϖne0 g hg0 hg1 hg n) ∧
      (∀ τ : K.absGal, (∀ y ∈ unramifiedClosure K, τ y = y) →
        ∀ lam ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n,
          e (τ lam) = τ (e lam)) ∧
      (∀ (a : 𝒪[K.carrier]) (lam : K.closure)
        (hlam : lam ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n)
        (hel : e lam ∈ iteratedLubinTateTorsionPoints K hq hϖmax hϖne0 g hg0 hg1 hg n),
        e (ltActPt K hq hπmax hπne0 f hf0 hf1 hf n lam hlam a)
          = ltActPt K hq hϖmax hϖne0 g hg0 hg1 hg n (e lam) hel a) := by
  obtain ⟨e, he, hbij⟩ := exists_bijOn_torsionPoints_of_dworkTheta K hq hπmax hπne0 f hf0 hf1 hf
    hϖmax hϖne0 g hg0 hg1 hg σ θ θ' hθ0 hθ'θ hint n
  refine ⟨e, hbij, ?_, ?_⟩
  · intro τ hτ lam hlam
    have hτlam : τ lam ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n :=
      algEquiv_mem_iteratedLubinTateTorsionPoints_of_mem K hq hπmax hπne0 f hf0 hf1 hf n τ lam
        hlam
    apply closureCompletionCoe_injective K
    rw [← he (τ lam) hτlam]
    exact coe_evalAt_algEquiv K τ hτ θ
      (norm_lt_one_of_mem_iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n lam hlam)
      (norm_lt_one_of_mem_iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n _ hτlam)
      (e lam) (he lam hlam)
  · intro a lam hlam hel
    have hμ := ltActPt_mem K hq hπmax hπne0 f hf0 hf1 hf n lam hlam a
    apply closureCompletionCoe_injective K
    rw [← he _ hμ]
    exact coe_evalAt_ltActPt K hq hπmax hπne0 f hf0 hf1 hf hϖmax hϖne0 g hg0 hg1 hg σ θ hθ0 hint
      n lam hlam
      (norm_lt_one_of_mem_iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n lam hlam)
      a
      (norm_lt_one_of_mem_iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n _ hμ)
      (e lam) hel
      (norm_lt_one_of_mem_iteratedLubinTateTorsionPoints K hq hϖmax hϖne0 g hg0 hg1 hg n _ hel)
      (he lam hlam)


/-! ## 8. ★★原始性の保存 -/

/-- ★★`𝒪_K`-線型な全単射は「原始的な捩れ点」を「原始的な捩れ点」へ移す。 -/
theorem transport_mem_psi (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal 𝒪[K.carrier]) 𝒪[K.carrier]]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField 𝒪[K.carrier]) pp]
    [Fintype (IsLocalRing.ResidueField 𝒪[K.carrier])]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField 𝒪[K.carrier]) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal 𝒪[K.carrier] = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries 𝒪[K.carrier]) (hf0 : PowerSeries.coeff 0 f = 0)
    (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue 𝒪[K.carrier]) f = PowerSeries.X ^ (pp ^ ff))
    {ϖ : 𝒪[K.carrier]} (hϖmax : IsLocalRing.maximalIdeal 𝒪[K.carrier] = Ideal.span {ϖ})
    (hϖne0 : ϖ ≠ 0)
    (g : PowerSeries 𝒪[K.carrier]) (hg0 : PowerSeries.coeff 0 g = 0)
    (hg1 : PowerSeries.coeff 1 g = ϖ)
    (hg : PowerSeries.map (IsLocalRing.residue 𝒪[K.carrier]) g = PowerSeries.X ^ (pp ^ ff))
    (n : ℕ) (hn : 1 ≤ n) (e : K.closure → K.closure)
    (hbij : Set.BijOn e ↑(iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n)
      ↑(iteratedLubinTateTorsionPoints K hq hϖmax hϖne0 g hg0 hg1 hg n))
    (hlin : ∀ (a : 𝒪[K.carrier]) (lam : K.closure)
      (hlam : lam ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n)
      (hel : e lam ∈ iteratedLubinTateTorsionPoints K hq hϖmax hϖne0 g hg0 hg1 hg n),
      e (ltActPt K hq hπmax hπne0 f hf0 hf1 hf n lam hlam a)
        = ltActPt K hq hϖmax hϖne0 g hg0 hg1 hg n (e lam) hel a)
    (lam : K.closure)
    (hlamψ : lam ∈ iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n hn) :
    e lam ∈ iteratedLubinTatePsiTorsionPoints K hq hϖmax hϖne0 g hg0 hg1 hg n hn := by
  rw [← iteratedLubinTateTorsionPoints_sdiff_eq_iteratedLubinTatePsiTorsionPoints
    K hq hπmax hπne0 f hf0 hf1 hf n hn, Finset.mem_sdiff] at hlamψ
  obtain ⟨hlamn, hlamn'⟩ := hlamψ
  have h0f : (0 : K.closure) ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n :=
    zero_mem_iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n
  have hel : e lam ∈ iteratedLubinTateTorsionPoints K hq hϖmax hϖne0 g hg0 hg1 hg n :=
    hbij.mapsTo hlamn
  have he0mem : e 0 ∈ iteratedLubinTateTorsionPoints K hq hϖmax hϖne0 g hg0 hg1 hg n :=
    hbij.mapsTo h0f
  have hz0 : ltActPt K hq hπmax hπne0 f hf0 hf1 hf n 0 h0f 0 = 0 := by
    show ((lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf n 0 h0f _ 0 :
      ↥(IntermediateField.adjoin K.carrier ({(0 : K.closure)} : Set K.closure))) : K.closure) = 0
    rw [lubinTateActionAtTorsionPoint_zero]
    rfl
  have hz0' : ltActPt K hq hϖmax hϖne0 g hg0 hg1 hg n (e 0) he0mem 0 = 0 := by
    show ((lubinTateActionAtTorsionPoint K hq hϖmax hϖne0 g hg0 hg1 hg n (e 0) he0mem _ 0 :
      ↥(IntermediateField.adjoin K.carrier ({e 0} : Set K.closure))) : K.closure) = 0
    rw [lubinTateActionAtTorsionPoint_zero]
    rfl
  have he0 : e 0 = 0 := by
    have := hlin 0 0 h0f he0mem
    rw [hz0, hz0'] at this
    exact this
  rw [← iteratedLubinTateTorsionPoints_sdiff_eq_iteratedLubinTatePsiTorsionPoints
    K hq hϖmax hϖne0 g hg0 hg1 hg n hn, Finset.mem_sdiff]
  refine ⟨hel, fun hcon => hlamn' ?_⟩
  have hpimem : π ^ (n - 1) ∈ IsLocalRing.maximalIdeal 𝒪[K.carrier] ^ (n - 1) := by
    rw [maximalIdeal_pow_eq_span_singleton_pi_pow K hπmax (n - 1)]
    exact Ideal.mem_span_singleton_self _
  have hgzero : ltActPt K hq hϖmax hϖne0 g hg0 hg1 hg n (e lam) hel (π ^ (n - 1)) = 0 := by
    have := (mem_iteratedLubinTateTorsionPoints_iff_forall_mem_maximalIdeal_pow
      K hq hϖmax hϖne0 g hg0 hg1 hg (e lam)
      (IntermediateField.mem_adjoin_simple_self K.carrier (e lam))
      (norm_lt_one_of_mem_iteratedLubinTateTorsionPoints K hq hϖmax hϖne0 g hg0 hg1 hg n _ hel)
      (n - 1)).mp hcon (π ^ (n - 1)) hpimem
    show ((lubinTateActionAtTorsionPoint K hq hϖmax hϖne0 g hg0 hg1 hg n (e lam) hel _
        (π ^ (n - 1)) : ↥(IntermediateField.adjoin K.carrier ({e lam} : Set K.closure)))
      : K.closure) = 0
    rw [← lubinTateActionAtPoint_eq_lubinTateActionAtTorsionPoint K hq hϖmax hϖne0 g hg0 hg1 hg
      (e lam) (IntermediateField.mem_adjoin_simple_self K.carrier (e lam))
      (norm_lt_one_of_mem_iteratedLubinTateTorsionPoints K hq hϖmax hϖne0 g hg0 hg1 hg n _ hel)
      n hel (π ^ (n - 1)), this]
    rfl
  have hfzero : ltActPt K hq hπmax hπne0 f hf0 hf1 hf n lam hlamn (π ^ (n - 1)) = 0 := by
    have hh := hlin (π ^ (n - 1)) lam hlamn hel
    rw [hgzero] at hh
    have := hbij.injOn (ltActPt_mem K hq hπmax hπne0 f hf0 hf1 hf n lam hlamn (π ^ (n - 1)))
      h0f (by rw [hh, he0])
    exact this
  refine (mem_iteratedLubinTateTorsionPoints_iff_lubinTateActionAtPoint_pi_pow_eq_zero
    K hq hπmax hπne0 f hf0 hf1 hf lam (IntermediateField.mem_adjoin_simple_self K.carrier lam)
    (norm_lt_one_of_mem_iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n lam hlamn)
    (n - 1)).mpr ?_
  rw [lubinTateActionAtPoint_eq_lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf
    lam (IntermediateField.mem_adjoin_simple_self K.carrier lam)
    (norm_lt_one_of_mem_iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n lam hlamn)
    n hlamn (π ^ (n - 1))]
  apply Subtype.ext
  apply Subtype.ext
  exact hfzero


/-! ## 9. ★★★★有限段の一致 `ρ_{f,n} = ρ_{g,n}`(慣性の上) -/

/-- ★★★★**有限段の一致** —— 輸送を持つとき、`ρ_{f,n}` と `ρ_{g,n}` は
慣性部分の上で一致する。 -/
theorem reciprocityMap_eq_of_transport (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal 𝒪[K.carrier]) 𝒪[K.carrier]]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField 𝒪[K.carrier]) pp]
    [Fintype (IsLocalRing.ResidueField 𝒪[K.carrier])]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField 𝒪[K.carrier]) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal 𝒪[K.carrier] = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries 𝒪[K.carrier]) (hf0 : PowerSeries.coeff 0 f = 0)
    (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue 𝒪[K.carrier]) f = PowerSeries.X ^ (pp ^ ff))
    {ϖ : 𝒪[K.carrier]} (hϖmax : IsLocalRing.maximalIdeal 𝒪[K.carrier] = Ideal.span {ϖ})
    (hϖne0 : ϖ ≠ 0)
    (g : PowerSeries 𝒪[K.carrier]) (hg0 : PowerSeries.coeff 0 g = 0)
    (hg1 : PowerSeries.coeff 1 g = ϖ)
    (hg : PowerSeries.map (IsLocalRing.residue 𝒪[K.carrier]) g = PowerSeries.X ^ (pp ^ ff))
    (n : ℕ) (hn : 1 ≤ n) (e : K.closure → K.closure)
    (hbij : Set.BijOn e ↑(iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n)
      ↑(iteratedLubinTateTorsionPoints K hq hϖmax hϖne0 g hg0 hg1 hg n))
    (heqv : ∀ τ : K.absGal, (∀ y ∈ unramifiedClosure K, τ y = y) →
      ∀ lam ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n,
        e (τ lam) = τ (e lam))
    (hlin : ∀ (a : 𝒪[K.carrier]) (lam : K.closure)
      (hlam : lam ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n)
      (hel : e lam ∈ iteratedLubinTateTorsionPoints K hq hϖmax hϖne0 g hg0 hg1 hg n),
      e (ltActPt K hq hπmax hπne0 f hf0 hf1 hf n lam hlam a)
        = ltActPt K hq hϖmax hϖne0 g hg0 hg1 hg n (e lam) hel a)
    (σ : K.absGal) (hσ : ∀ y ∈ unramifiedClosure K, σ y = y) (u : (𝒪[K.carrier])ˣ)
    (x : K.closure)
    (hxψ : x ∈ iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n hn)
    (hxn : x ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n)
    (hmem : x ∈ IntermediateField.adjoin K.carrier ({x} : Set K.closure))
    (hrec : reciprocityMap K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem σ
      = QuotientGroup.mk u)
    (y : K.closure)
    (hyψ : y ∈ iteratedLubinTatePsiTorsionPoints K hq hϖmax hϖne0 g hg0 hg1 hg n hn)
    (hyn : y ∈ iteratedLubinTateTorsionPoints K hq hϖmax hϖne0 g hg0 hg1 hg n)
    (hymem : y ∈ IntermediateField.adjoin K.carrier ({y} : Set K.closure)) :
    reciprocityMap K hq hϖmax hϖne0 g hg0 hg1 hg n hn y hyψ hyn hymem σ
      = QuotientGroup.mk u := by
  have hel : e x ∈ iteratedLubinTateTorsionPoints K hq hϖmax hϖne0 g hg0 hg1 hg n :=
    hbij.mapsTo hxn
  have helψ : e x ∈ iteratedLubinTatePsiTorsionPoints K hq hϖmax hϖne0 g hg0 hg1 hg n hn :=
    transport_mem_psi K hq hπmax hπne0 f hf0 hf1 hf hϖmax hϖne0 g hg0 hg1 hg n hn e hbij hlin
      x hxψ
  have hspec := reciprocityMap_spec K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem σ
  rw [hrec, unitActionQuotientLift_mk] at hspec
  have hx' : σ x = ltActPt K hq hπmax hπne0 f hf0 hf1 hf n x hxn (u : 𝒪[K.carrier]) := hspec.symm
  have h3 : σ (e x)
      = ltActPt K hq hϖmax hϖne0 g hg0 hg1 hg n (e x) hel (u : 𝒪[K.carrier]) := by
    rw [← hlin (u : 𝒪[K.carrier]) x hxn hel, ← hx']
    exact (heqv σ hσ x hxn).symm
  have hmk : reciprocityMap K hq hϖmax hϖne0 g hg0 hg1 hg n hn (e x) helψ hel
      (IntermediateField.mem_adjoin_simple_self K.carrier (e x)) σ = QuotientGroup.mk u :=
    reciprocityMap_eq_mk_of_apply_eq K hq hϖmax hϖne0 g hg0 hg1 hg n hn (e x) helψ hel
      (IntermediateField.mem_adjoin_simple_self K.carrier (e x)) u σ h3
  rw [reciprocityMap_point_indep K hq hϖmax hϖne0 g hg0 hg1 hg n hn y (e x) hyψ hyn hymem
    helψ hel (IntermediateField.mem_adjoin_simple_self K.carrier (e x)) σ]
  exact hmk


/-! ## 10. ★★★★★極限への持ち上げ -/

/-- ★★★★★**極限への持ち上げ** —— 各段の輸送があれば `ρ_g(σ) = ρ_f(σ)`。 -/
theorem reciprocityUnits_eq_of_transport (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal 𝒪[K.carrier]) 𝒪[K.carrier]]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField 𝒪[K.carrier]) pp]
    [Fintype (IsLocalRing.ResidueField 𝒪[K.carrier])]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField 𝒪[K.carrier]) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal 𝒪[K.carrier] = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries 𝒪[K.carrier]) (hf0 : PowerSeries.coeff 0 f = 0)
    (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue 𝒪[K.carrier]) f = PowerSeries.X ^ (pp ^ ff))
    {ϖ : 𝒪[K.carrier]} (hϖmax : IsLocalRing.maximalIdeal 𝒪[K.carrier] = Ideal.span {ϖ})
    (hϖne0 : ϖ ≠ 0)
    (g : PowerSeries 𝒪[K.carrier]) (hg0 : PowerSeries.coeff 0 g = 0)
    (hg1 : PowerSeries.coeff 1 g = ϖ)
    (hg : PowerSeries.map (IsLocalRing.residue 𝒪[K.carrier]) g = PowerSeries.X ^ (pp ^ ff))
    (σ : K.absGal) (hσ : ∀ y ∈ unramifiedClosure K, σ y = y)
    (htrans : ∀ n : ℕ, ∃ e : K.closure → K.closure,
      Set.BijOn e ↑(iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n)
        ↑(iteratedLubinTateTorsionPoints K hq hϖmax hϖne0 g hg0 hg1 hg n) ∧
      (∀ τ : K.absGal, (∀ y ∈ unramifiedClosure K, τ y = y) →
        ∀ lam ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n,
          e (τ lam) = τ (e lam)) ∧
      (∀ (a : 𝒪[K.carrier]) (lam : K.closure)
        (hlam : lam ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n)
        (hel : e lam ∈ iteratedLubinTateTorsionPoints K hq hϖmax hϖne0 g hg0 hg1 hg n),
        e (ltActPt K hq hπmax hπne0 f hf0 hf1 hf n lam hlam a)
          = ltActPt K hq hϖmax hϖne0 g hg0 hg1 hg n (e lam) hel a)) :
    reciprocityUnits K hq hϖmax hϖne0 g hg0 hg1 hg σ
      = reciprocityUnits K hq hπmax hπne0 f hf0 hf1 hf σ := by
  have hspan : ∀ m : ℕ, Ideal.span ({ϖ ^ m} : Set 𝒪[K.carrier])
      = Ideal.span ({π ^ m} : Set 𝒪[K.carrier]) := by
    intro m
    rw [← maximalIdeal_pow_eq_span_singleton_pi_pow K hϖmax m,
      maximalIdeal_pow_eq_span_singleton_pi_pow K hπmax m]
  apply (unitsEquivCompatibleUnits K hϖmax).injective
  apply Subtype.ext
  funext m
  show unitReductionQuotientMap K ϖ m _ = unitReductionQuotientMap K ϖ m _
  cases m with
  | zero => exact @Subsingleton.elim _ (subsingleton_units_quotient_pow_zero K _) _ _
  | succ k =>
    obtain ⟨u, hu⟩ := QuotientGroup.mk_surjective
      (reciprocityMap K hq hπmax hπne0 f hf0 hf1 hf (k + 1) (by omega)
        (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf k).pt
        (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf k).hψ
        (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf k).hn
        (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf k).hmem σ)
    obtain ⟨e, hbij, heqv, hlin⟩ := htrans (k + 1)
    have hgk := reciprocityMap_eq_of_transport K hq hπmax hπne0 f hf0 hf1 hf
      hϖmax hϖne0 g hg0 hg1 hg (k + 1) (by omega) e hbij heqv hlin σ hσ u
      (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf k).pt
      (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf k).hψ
      (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf k).hn
      (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf k).hmem hu.symm
      (psiGenSeq K hq hϖmax hϖne0 g hg0 hg1 hg k).pt
      (psiGenSeq K hq hϖmax hϖne0 g hg0 hg1 hg k).hψ
      (psiGenSeq K hq hϖmax hϖne0 g hg0 hg1 hg k).hn
      (psiGenSeq K hq hϖmax hϖne0 g hg0 hg1 hg k).hmem
    have hsub : (u : 𝒪[K.carrier])
        - ((reciprocityUnits K hq hπmax hπne0 f hf0 hf1 hf σ : (𝒪[K.carrier])ˣ) :
            𝒪[K.carrier]) ∈ Ideal.span ({π ^ (k + 1)} : Set 𝒪[K.carrier]) := by
      rw [← Ideal.Quotient.mk_eq_mk_iff_sub_mem]
      have h1 : unitReductionQuotientMap K π (k + 1) u
          = unitReductionQuotientMap K π (k + 1)
              (reciprocityUnits K hq hπmax hπne0 f hf0 hf1 hf σ) := by
        rw [unitReductionQuotientMap_reciprocityUnits K hq hπmax hπne0 f hf0 hf1 hf σ (k + 1),
          ← principalUnitsQuotientEquiv_apply_mk K hπmax (k + 1) (by omega) u, hu]
        rfl
      exact congrArg Units.val h1
    refine Eq.trans (unitReductionQuotientMap_reciprocityUnits K hq hϖmax hϖne0 g hg0 hg1 hg
      σ (k + 1)) ?_
    refine Eq.trans (congrArg (principalUnitsQuotientEquiv K hϖmax (k + 1) (by omega)) hgk) ?_
    refine Eq.trans (principalUnitsQuotientEquiv_apply_mk K hϖmax (k + 1) (by omega) u) ?_
    apply Units.ext
    show Ideal.Quotient.mk (Ideal.span ({ϖ ^ (k + 1)} : Set 𝒪[K.carrier])) (u : 𝒪[K.carrier])
      = Ideal.Quotient.mk (Ideal.span ({ϖ ^ (k + 1)} : Set 𝒪[K.carrier]))
        ((reciprocityUnits K hq hπmax hπne0 f hf0 hf1 hf σ : (𝒪[K.carrier])ˣ) : 𝒪[K.carrier])
    rw [Ideal.Quotient.mk_eq_mk_iff_sub_mem, hspan (k + 1)]
    exact hsub


/-! ## 11. ★★★★★★★★穴の充填と pGC Proposition 1.1

原文 (pGC p.3):
> The cyclotomic character χ : Γ_K → Z[bb]_p^× can be recovered entirely
> group-theoretically from Γ_K.

★★**ここで `ReciprocityDatumIndependenceOnTorsion` が無条件で出る。** -/

/-- ★★★★★★**慣性の上での素元非依存性**(原典 Corollary 4.9 の慣性部分)。 -/
theorem reciprocityUnits_eq_of_isUnit_mul (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal 𝒪[K.carrier]) 𝒪[K.carrier]]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField 𝒪[K.carrier]) pp]
    [Fintype (IsLocalRing.ResidueField 𝒪[K.carrier])]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField 𝒪[K.carrier]) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal 𝒪[K.carrier] = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries 𝒪[K.carrier]) (hf0 : PowerSeries.coeff 0 f = 0)
    (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue 𝒪[K.carrier]) f = PowerSeries.X ^ (pp ^ ff))
    {u : 𝒪[K.carrier]} (hu : IsUnit u)
    (hϖmax : IsLocalRing.maximalIdeal 𝒪[K.carrier] = Ideal.span {u * π}) (hϖne0 : u * π ≠ 0)
    (g : PowerSeries 𝒪[K.carrier]) (hg0 : PowerSeries.coeff 0 g = 0)
    (hg1 : PowerSeries.coeff 1 g = u * π)
    (hg : PowerSeries.map (IsLocalRing.residue 𝒪[K.carrier]) g = PowerSeries.X ^ (pp ^ ff))
    (σ : K.absGal) (hσ : ∀ y ∈ unramifiedClosure K, σ y = y) :
    reciprocityUnits K hq hϖmax hϖne0 g hg0 hg1 hg σ
      = reciprocityUnits K hq hπmax hπne0 f hf0 hf1 hf σ := by
  obtain ⟨τ, -, -, -, hstep2⟩ := exists_arithFrobenius_isCoherent_dworkThetaStep2 K hq
  obtain ⟨θ, θ', hθ0, -, -, hθ'θ, -, hint⟩ := hstep2 π hπmax f hf0 hf1 hf u hu g hg0 hg1 hg
  refine reciprocityUnits_eq_of_transport K hq hπmax hπne0 f hf0 hf1 hf hϖmax hϖne0 g hg0 hg1 hg
    σ hσ (fun n => ?_)
  exact exists_transport_of_dworkTheta K hq hπmax hπne0 f hf0 hf1 hf hϖmax hϖne0 g hg0 hg1 hg
    τ θ θ' ((PowerSeries.coeff_zero_eq_constantCoeff_apply θ).symm.trans hθ0) hθ'θ hint n

/-- ★★★★★★★★**残っていた穴が埋まった**。 -/
theorem reciprocityDatumIndependenceOnTorsion_holds :
    ReciprocityDatumIndependenceOnTorsion p := by
  intro K _ pp _ _ ff hq π hπmax hπne0 f hf0 hf1 hf π' hπmax' hπne0' f' hf0' hf1' hf' m θ htor
  have hassoc : Associated π π' := by
    rw [← Ideal.span_singleton_eq_span_singleton (α := 𝒪[K.carrier])]
    rw [← hπmax, hπmax']
  obtain ⟨v, hv⟩ := hassoc
  have hπ'eq : π' = (v : 𝒪[K.carrier]) * π := by rw [← hv, mul_comm]
  subst hπ'eq
  exact congrArg (fun w : (𝒪[K.carrier])ˣ => (w : 𝒪[K.carrier]))
    (reciprocityUnits_eq_of_isUnit_mul K hq hπmax hπne0 f hf0 hf1 hf v.isUnit
      hπmax' hπne0' f' hf0' hf1' hf' θ
      (fun _ hy => apply_eq_self_of_mem_unramifiedClosure_of_torsion K m θ htor hy))

/-- ★★★★★★★★★★**壁が無条件で出た** —— `ArtinUnitEquivariance p`。 -/
theorem artinUnitEquivariance_holds : ArtinUnitEquivariance p :=
  artinUnitEquivariance_of_reciprocityDatumIndependenceOnTorsion
    reciprocityDatumIndependenceOnTorsion_holds

/-- ★★★★★★★★★★★★**[pGC] Proposition 1.1(無条件)**。 -/
theorem cyclotomicCharacter_recoverable_holds :
    (cyclotomicCharacterObject (p := p)).RecoverableFromAbsGal :=
  cyclotomicCharacter_recoverable_of_artinUnitEquivariance artinUnitEquivariance_holds

end ABC3.Found.PGC
