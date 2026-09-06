import ABC3.Found.PGC.UnramifiedZhat
import Mathlib.Analysis.Normed.Algebra.Ultra
import Mathlib.Analysis.Normed.Module.Completion
import Mathlib.Topology.Algebra.Valued.NormedValued
import Mathlib.Topology.Algebra.Valued.ValuedField

/-!
# `K̂^{ur}`(不分岐閉包の完備化)と、その上の Frobenius(`sorry` 無し)

経路 Λ の節点 Λ5 の後半。Λ5 の前半(`Found/PGC/UnramifiedZhat.lean`)は
`Gal(K^ur/K) ≃ₜ* Ẑ` を示した。本ファイルは **Λ6(`Art_π` の π 非依存性
= Dwork の補題)の舞台**を建てる:

```
K̂^{ur} := (K^ur の完備化)   𝒪_{K̂^{ur}}   Frobenius Φ : K̂^{ur} ≃ₐ[K] K̂^{ur}
```

## mathlib の在庫でどこまで済んだか(2026-09-06 の実測)

★**見積 800–1,500 行に対し実際は約 330 行**。理由は、必要な構造がすべて
mathlib の既製品に乗ったこと。Λ1・Λ3・Λ4 と同じ現象である。

| 要るもの | mathlib の在庫 |
|---|---|
| `K^ur` のノルム | `spectralNorm.normedField`(`K.closure` 経由、`LocalFieldNorm.lean` が既に置いている) |
| 完備化 | `UniformSpace.Completion` + `Valued.completable`(付値体は completable) |
| 完備化の体構造 | `UniformSpace.Completion.instField`(`CompletableTopField` から) |
| 完備化のノルム | `[NormedField A] [CompletableTopField A] : NormedField (Completion A)` |
| 完備化の付値 | `NormedField.toValued`(超距離ノルム → `Valued _ ℝ≥0`) |
| 整数環 | `Valuation.valuationSubring` / `Valuation.integer.integers` |
| 自己同型の延長 | `UniformSpace.Completion.mapRingEquiv` |
| Galois が等長 | `spectralNorm_eq_of_equiv`(`σ` はスペクトルノルムを保つ) |

雛形は `Mathlib/NumberTheory/Padics/Complex.lean`(`ℂ_[p]` = `PadicAlgCl p` の
完備化)。本ファイルは同じ骨組みを「代数閉包」でなく「不分岐閉包」で回す。

## 到達点

1. **`K̂^{ur}` の構成**(`unramifiedCompletion`)。完備なノルム体・付値体で、
   `K` 上のノルム代数。`K^ur` は稠密(`denseRange_coe_unramifiedCompletion`)。
2. **`𝒪_{K̂^{ur}}`**(`unramifiedCompletionInt`)。`ValuationSubring` なので
   `IsLocalRing` / `ValuationRing` が付く。極大イデアルは `‖z‖ < 1`
   (`mem_maximalIdeal_unramifiedCompletionInt`)。
   さらに **`𝓀_{K^{ur}} ≅ 𝓀_{K̂^{ur}}`**(`residueFieldEquiv`)——完備化は
   剰余体を変えない。
3. **Frobenius の持ち上げ**(`exists_frobenius_lift`)。`Gal(K^ur/K)` の元は
   スペクトルノルムに関して等長なので一意に連続延長し、
   `Gal(K^ur/K) →* Aut(K̂^{ur}/K)` は**単射**
   (`unramGalCompletionHom_injective`)。`UnramifiedZhat.lean` の整合的
   Frobenius を延長したものが `frobeniusCompletion` で、その
   `𝒪_{K̂^{ur}}` への制限が `frobeniusCompletionInt`(Dwork の補題の舞台)。

## 到達していない点(正直な記録)

* **Dwork の補題**(`Φ(ξ)/ξ = u` を解く)は**入っていない**。逐次近似の
  材料(完備性・`exists_norm_le_one_norm_sub_lt`・`frobeniusCompletionInt`)は
  揃ったが、各段の補正を作るのに
  「剰余体 `𝓀_{K̂^{ur}}` が代数閉体 `𝔽̄_q` で、Frobenius がその上で `x ↦ x^q`」
  が要る。**そこは本ファイルにも `UnramifiedZhat.lean` にも無い**。
* すなわち Λ6 の手前に**新しい節点**が要る:
  「`𝓀_{K^{ur}} ≅ 𝔽̄_q` かつ整合的 Frobenius の剰余体上の作用が `x ↦ x^q`」。
  `UnramifiedZhat.lean` は意図的に剰余体へ降りずに `Ẑ` を作ったので
  (同ファイルの逸脱記録)、この同定は**未証明のまま残っている**。
  `residueFieldEquiv` はその節点の半分(完備化側と `K^ur` 側の一致)を
  先に片付けたものである。

## ★設計上の注意(守ったこと)

* **`Valued` を `K^ur` の instance にしていない**。もし instance にすると
  `Valued.valuedCompletion` が `K̂^{ur}` に**もう一つの** `Valued` を作り、
  `NormedField.toValued` と菱形になる。`K^ur` 側では
  `NormedField.valuation`(束ねない `Valuation`)だけを使い、
  `CompletableTopField` は `letI` の中で `Valued.completable` を借りる。
  この選び方のおかげで `Valued.v z = ‖z‖₊` が `K̂^{ur}` 上で `rfl` になる。
* **`inertia` を経由していない**。不分岐側は `unramifiedClosure K` という体として
  直接扱っており、`Interface` の `SubgroupCorrespondence` / `ResidueCardinality` は
  本ファイルの主張にも証明にも現れない(Corollary 1.3 との循環を避けるため)。
* **`Abelianization` を使っていない**。
* **結論に自由なパラメータを出していない**——`exists_frobenius_lift` の型は
  `K` にしか依存せず、生成元 `σ` は `∃` の内側にある。
  `frobeniusCompletion` の `σ` は `Classical.choose` に閉じ込めてある。
* **`sorry` 本体の `def` を作っていない**(D21)。定義できないもの
  (剰余体の Frobenius 同定)は**書かずに上の「到達していない点」に記録した**。

## 逸脱(記録)

* 古典的には `K̂^{ur}` は「離散付値環 `𝒪_{K^{ur}}` の `𝔪`-進完備化の商体」として
  作ることが多い。本ファイルは**ノルム(スペクトルノルム)の完備化**として作る。
  `K^ur/K` は代数拡大なので値群は `ℚ` の離散部分群に入り、両者は同じものだが、
  一致の証明は付けていない(下流で必要になった時点で節点を立てること)。
* 「Frobenius」は本ファイルでも `UnramifiedZhat.lean` の
  `exists_coherentFrobenius` が返す元を指す。剰余体上の `x ↦ x^q` との一致は
  **主張していない**。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC
open scoped NNReal

variable {p : ℕ} [Fact p.Prime]

/-! ## 1. `K^ur` のノルムと、Galois 作用が等長であること -/

/-- `K^ur` のノルム(`K.closure` から誘導したもの)は `K` 上のスペクトルノルム。 -/
theorem norm_eq_spectralNorm_unramifiedClosure (K : PAdicLocalField p)
    (x : ↥(unramifiedClosure K)) : ‖x‖ = spectralNorm K.carrier ↥(unramifiedClosure K) x :=
  (spectralNorm.eq_of_tower (K := K.carrier) (L := K.closure)
    (E := ↥(unramifiedClosure K)) x).symm

/-- **`Gal(K^ur/K)` の元は等長**。スペクトルノルムは最小多項式だけで決まり、
`σ` は最小多項式を変えないから。 -/
theorem norm_unramGal (K : PAdicLocalField p) (σ : unramGal K) (x : ↥(unramifiedClosure K)) :
    ‖σ x‖ = ‖x‖ := by
  rw [norm_eq_spectralNorm_unramifiedClosure, norm_eq_spectralNorm_unramifiedClosure,
    ← spectralNorm_eq_of_equiv σ x]

theorem isometry_unramGal (K : PAdicLocalField p) (σ : unramGal K) :
    Isometry (σ : ↥(unramifiedClosure K) → ↥(unramifiedClosure K)) :=
  Isometry.of_dist_eq fun x y => by
    rw [dist_eq_norm, dist_eq_norm, ← map_sub, norm_unramGal]

/-- `K^ur` は completable(付値体だから)。

★`Valued` を instance にせず `letI` で借りているのは、`Valued.valuedCompletion` が
`K̂^{ur}` に第二の `Valued` を作るのを防ぐため(モジュール docstring 参照)。 -/
noncomputable scoped instance completableTopField_unramifiedClosure (K : PAdicLocalField p) :
    CompletableTopField ↥(unramifiedClosure K) :=
  letI : Valued ↥(unramifiedClosure K) NNReal := NormedField.toValued
  Valued.completable

/-! ## 2. `K̂^{ur}` -/

/-- **`K̂^{ur}`** —— 不分岐閉包 `K^ur` のスペクトルノルムによる完備化。 -/
abbrev unramifiedCompletion (K : PAdicLocalField p) : Type :=
  UniformSpace.Completion ↥(unramifiedClosure K)

theorem norm_algebraMap_unramifiedClosure (K : PAdicLocalField p) (x : K.carrier) :
    ‖algebraMap K.carrier ↥(unramifiedClosure K) x‖ = ‖x‖ := by
  rw [norm_eq_spectralNorm_unramifiedClosure]
  exact spectralNorm_extends x

noncomputable scoped instance isUltrametricDist_unramifiedCompletion (K : PAdicLocalField p) :
    IsUltrametricDist (unramifiedCompletion K) := IsUltrametricDist.of_normedAlgebra K.carrier

/-- `K̂^{ur}` の付値。**`Valued.v z = ‖z‖₊` が `rfl`** になるように、
完備化の `Valued.valuedCompletion` ではなくノルムから作る。 -/
noncomputable scoped instance valued_unramifiedCompletion (K : PAdicLocalField p) :
    Valued (unramifiedCompletion K) NNReal := NormedField.toValued

theorem norm_coe_unramifiedCompletion (K : PAdicLocalField p) (x : ↥(unramifiedClosure K)) :
    ‖(x : unramifiedCompletion K)‖ = ‖x‖ := UniformSpace.Completion.norm_coe x

theorem norm_algebraMap_unramifiedCompletion (K : PAdicLocalField p) (x : K.carrier) :
    ‖algebraMap K.carrier (unramifiedCompletion K) x‖ = ‖x‖ := by
  rw [UniformSpace.Completion.algebraMap_def, norm_coe_unramifiedCompletion,
    norm_algebraMap_unramifiedClosure]

theorem valuation_unramifiedCompletion (K : PAdicLocalField p) (z : unramifiedCompletion K) :
    Valued.v z = ‖z‖₊ := rfl

theorem valuation_coe_unramifiedCompletion (K : PAdicLocalField p) (x : ↥(unramifiedClosure K)) :
    Valued.v (x : unramifiedCompletion K) = ‖x‖₊ := by
  rw [valuation_unramifiedCompletion]
  ext
  exact norm_coe_unramifiedCompletion K x

noncomputable scoped instance nontriviallyNormedField_unramifiedCompletion
    (K : PAdicLocalField p) : NontriviallyNormedField (unramifiedCompletion K) where
  non_trivial := by
    obtain ⟨x, hx⟩ := NontriviallyNormedField.non_trivial (α := K.carrier)
    exact ⟨algebraMap K.carrier (unramifiedCompletion K) x,
      by rwa [norm_algebraMap_unramifiedCompletion]⟩

noncomputable scoped instance charZero_unramifiedCompletion (K : PAdicLocalField p) :
    CharZero (unramifiedCompletion K) := by
  haveI : CharZero K.carrier :=
    charZero_of_injective_algebraMap (algebraMap ℚ_[p] K.carrier).injective
  exact (RingHom.charZero_iff
    (algebraMap K.carrier (unramifiedCompletion K)).injective).mp inferInstance

/-- `K^ur` は `K̂^{ur}` で稠密。 -/
theorem denseRange_coe_unramifiedCompletion (K : PAdicLocalField p) :
    DenseRange (fun x : ↥(unramifiedClosure K) => (x : unramifiedCompletion K)) :=
  UniformSpace.Completion.denseRange_coe

/-- **任意精度の近似**。 -/
theorem exists_norm_sub_lt (K : PAdicLocalField p) (z : unramifiedCompletion K) {ε : ℝ}
    (hε : 0 < ε) : ∃ x : ↥(unramifiedClosure K), ‖z - (x : unramifiedCompletion K)‖ < ε := by
  obtain ⟨-, ⟨x, rfl⟩, hx⟩ :=
    Metric.mem_closure_iff.mp
      ((denseRange_coe_unramifiedCompletion K).closure_eq ▸ Set.mem_univ z) ε hε
  exact ⟨x, by rwa [← dist_eq_norm]⟩

/-- **整数の近似**——`‖z‖ ≤ 1` なら近似元も `𝒪_{K^{ur}}` の中に取れる。

★Λ6(Dwork の補題)の逐次近似はこの形を使う。

退化の自己検査:`hz` を落とすと `‖x‖ ≤ 1` は**偽**——`‖z‖ > 1` なら
`ε` を小さく取ったとき `‖x‖ = ‖z‖ > 1` になる。 -/
theorem exists_norm_le_one_norm_sub_lt (K : PAdicLocalField p) (z : unramifiedCompletion K)
    (hz : ‖z‖ ≤ 1) {ε : ℝ} (hε : 0 < ε) :
    ∃ x : ↥(unramifiedClosure K), ‖x‖ ≤ 1 ∧ ‖z - (x : unramifiedCompletion K)‖ < ε := by
  obtain ⟨x, hx⟩ := exists_norm_sub_lt K z (lt_min hε one_pos)
  refine ⟨x, ?_, lt_of_lt_of_le hx (min_le_left _ _)⟩
  have h1 : ‖(x : unramifiedCompletion K)‖ ≤ 1 := by
    have hmax := IsUltrametricDist.norm_add_le_max z (-(z - (x : unramifiedCompletion K)))
    simp only [neg_sub, add_sub_cancel] at hmax
    refine hmax.trans (max_le hz ?_)
    rw [norm_sub_rev]
    exact le_of_lt (lt_of_lt_of_le hx (min_le_right _ _))
  rwa [norm_coe_unramifiedCompletion] at h1

/-! ## 3. `𝒪_{K̂^{ur}}` -/

/-- **`𝒪_{K̂^{ur}}`** —— `K̂^{ur}` の整数環(付値 `≤ 1` の元)。 -/
noncomputable def unramifiedCompletionInt (K : PAdicLocalField p) :
    ValuationSubring (unramifiedCompletion K) :=
  (Valued.v : Valuation (unramifiedCompletion K) NNReal).valuationSubring

@[simp] theorem mem_unramifiedCompletionInt (K : PAdicLocalField p) (z : unramifiedCompletion K) :
    z ∈ unramifiedCompletionInt K ↔ ‖z‖ ≤ 1 := by
  rw [unramifiedCompletionInt, Valuation.mem_valuationSubring_iff]
  exact ⟨fun h => by exact_mod_cast h, fun h => by exact_mod_cast h⟩

theorem unramifiedCompletionInt_integers (K : PAdicLocalField p) :
    (Valued.v : Valuation (unramifiedCompletion K) NNReal).Integers
      (unramifiedCompletionInt K) :=
  Valuation.integer.integers _

/-! ## 4. `Gal(K^ur/K)` の `K̂^{ur}` への延長 -/

/-- `σ ∈ Gal(K^ur/K)` の完備化への延長(環同型として)。 -/
noncomputable def unramGalCompletionRingEquiv (K : PAdicLocalField p) (σ : unramGal K) :
    unramifiedCompletion K ≃+* unramifiedCompletion K :=
  UniformSpace.Completion.mapRingEquiv σ.toRingEquiv (isometry_unramGal K σ).continuous
    (isometry_unramGal K σ.symm).continuous

theorem unramGalCompletionRingEquiv_coe (K : PAdicLocalField p) (σ : unramGal K)
    (x : ↥(unramifiedClosure K)) :
    unramGalCompletionRingEquiv K σ (x : unramifiedCompletion K)
      = ((σ x : ↥(unramifiedClosure K)) : unramifiedCompletion K) :=
  UniformSpace.Completion.map_coe (isometry_unramGal K σ).uniformContinuous x

theorem continuous_unramGalCompletionRingEquiv (K : PAdicLocalField p) (σ : unramGal K) :
    Continuous (unramGalCompletionRingEquiv K σ) := UniformSpace.Completion.continuous_map

/-- **`σ ∈ Gal(K^ur/K)` の `K̂^{ur}` への延長**(`K`-代数同型)。 -/
noncomputable def unramGalCompletion (K : PAdicLocalField p) (σ : unramGal K) :
    unramifiedCompletion K ≃ₐ[K.carrier] unramifiedCompletion K :=
  AlgEquiv.ofRingEquiv (f := unramGalCompletionRingEquiv K σ) fun r => by
    rw [UniformSpace.Completion.algebraMap_def, unramGalCompletionRingEquiv_coe]
    congr 1
    exact σ.commutes r

@[simp] theorem unramGalCompletion_coe (K : PAdicLocalField p) (σ : unramGal K)
    (x : ↥(unramifiedClosure K)) :
    unramGalCompletion K σ (x : unramifiedCompletion K)
      = ((σ x : ↥(unramifiedClosure K)) : unramifiedCompletion K) :=
  unramGalCompletionRingEquiv_coe K σ x

theorem continuous_unramGalCompletion (K : PAdicLocalField p) (σ : unramGal K) :
    Continuous (unramGalCompletion K σ) := continuous_unramGalCompletionRingEquiv K σ

/-- 連続な 2 つの `K`-代数同型は `K^ur` 上で一致すれば等しい(稠密性)。 -/
theorem unramGalCompletion_ext (K : PAdicLocalField p)
    {f g : unramifiedCompletion K ≃ₐ[K.carrier] unramifiedCompletion K}
    (hf : Continuous f) (hg : Continuous g)
    (h : ∀ x : ↥(unramifiedClosure K),
      f (x : unramifiedCompletion K) = g (x : unramifiedCompletion K)) : f = g :=
  AlgEquiv.ext (congrFun (UniformSpace.Completion.ext hf hg h))

/-- **★★★★★★★★★★(Λ5)`Gal(K^ur/K) → Aut(K̂^{ur}/K)`**。 -/
noncomputable def unramGalCompletionHom (K : PAdicLocalField p) :
    unramGal K →* (unramifiedCompletion K ≃ₐ[K.carrier] unramifiedCompletion K) where
  toFun := unramGalCompletion K
  map_one' := unramGalCompletion_ext K (continuous_unramGalCompletion K 1) continuous_id (by simp)
  map_mul' σ τ := unramGalCompletion_ext K (continuous_unramGalCompletion K (σ * τ))
    ((continuous_unramGalCompletion K σ).comp (continuous_unramGalCompletion K τ)) (by simp)

@[simp] theorem unramGalCompletionHom_apply (K : PAdicLocalField p) (σ : unramGal K) :
    unramGalCompletionHom K σ = unramGalCompletion K σ := rfl

/-- **延長は単射**——`K^ur ↪ K̂^{ur}` が単射だから。

退化の自己検査:全射ではない。`K̂^{ur}` の連続でない `K`-代数自己同型が
あるかどうかは本ファイルの射程外で、`unramGalCompletionHom` の像は
少なくとも連続なものの中に入る。 -/
theorem unramGalCompletionHom_injective (K : PAdicLocalField p) :
    Function.Injective (unramGalCompletionHom K) := by
  intro σ τ h
  ext x
  have hx : ((σ x : ↥(unramifiedClosure K)) : unramifiedCompletion K)
      = ((τ x : ↥(unramifiedClosure K)) : unramifiedCompletion K) := by
    rw [← unramGalCompletion_coe, ← unramGalCompletion_coe, ← unramGalCompletionHom_apply,
      ← unramGalCompletionHom_apply, h]
  exact congrArg Subtype.val
    (UniformSpace.Completion.coe_injective ↥(unramifiedClosure K) hx)

/-- **延長も等長**。 -/
theorem norm_unramGalCompletion (K : PAdicLocalField p) (σ : unramGal K)
    (z : unramifiedCompletion K) : ‖unramGalCompletion K σ z‖ = ‖z‖ := by
  have h : (fun w : unramifiedCompletion K => ‖unramGalCompletion K σ w‖)
      = (fun w : unramifiedCompletion K => ‖w‖) :=
    UniformSpace.Completion.ext
      (continuous_norm.comp (continuous_unramGalCompletion K σ)) continuous_norm
      (fun x => by rw [unramGalCompletion_coe, norm_coe_unramifiedCompletion,
        norm_coe_unramifiedCompletion, norm_unramGal])
  exact congrFun h z

theorem isometry_unramGalCompletion (K : PAdicLocalField p) (σ : unramGal K) :
    Isometry (unramGalCompletion K σ) :=
  Isometry.of_dist_eq fun x y => by
    rw [dist_eq_norm, dist_eq_norm, ← map_sub, norm_unramGalCompletion]

/-- **延長は付値を保つ**。 -/
theorem valuation_unramGalCompletion (K : PAdicLocalField p) (σ : unramGal K)
    (z : unramifiedCompletion K) : Valued.v (unramGalCompletion K σ z) = Valued.v z := by
  rw [valuation_unramifiedCompletion, valuation_unramifiedCompletion]
  ext
  exact norm_unramGalCompletion K σ z

/-- **`𝒪_{K̂^{ur}}` は Galois 作用で保たれる**。 -/
@[simp] theorem unramGalCompletion_mem_int (K : PAdicLocalField p) (σ : unramGal K)
    (z : unramifiedCompletion K) :
    unramGalCompletion K σ z ∈ unramifiedCompletionInt K ↔ z ∈ unramifiedCompletionInt K := by
  rw [mem_unramifiedCompletionInt, mem_unramifiedCompletionInt, norm_unramGalCompletion]

/-! ## 5. Frobenius -/

/-- 整合的 Frobenius(`UnramifiedZhat.lean` の `exists_coherentFrobenius` を選択したもの)。 -/
noncomputable def coherentFrobenius (K : PAdicLocalField p) : unramGal K :=
  (exists_coherentFrobenius K).choose

theorem coherentFrobenius_spec (K : PAdicLocalField p) (N : ℕ) (hN : N ≠ 0) :
    coherentFrobenius K ∈ unramLevelGeneratorSet K N :=
  (exists_coherentFrobenius K).choose_spec N hN

/-- **`K̂^{ur}` の Frobenius**——`Gal(K^ur/K)` の位相的生成元を完備化へ延長したもの。 -/
noncomputable def frobeniusCompletion (K : PAdicLocalField p) :
    unramifiedCompletion K ≃ₐ[K.carrier] unramifiedCompletion K :=
  unramGalCompletion K (coherentFrobenius K)

@[simp] theorem frobeniusCompletion_coe (K : PAdicLocalField p) (x : ↥(unramifiedClosure K)) :
    frobeniusCompletion K (x : unramifiedCompletion K)
      = ((coherentFrobenius K x : ↥(unramifiedClosure K)) : unramifiedCompletion K) :=
  unramGalCompletion_coe K _ x

theorem continuous_frobeniusCompletion (K : PAdicLocalField p) :
    Continuous (frobeniusCompletion K) := continuous_unramGalCompletion K _

theorem norm_frobeniusCompletion (K : PAdicLocalField p) (z : unramifiedCompletion K) :
    ‖frobeniusCompletion K z‖ = ‖z‖ := norm_unramGalCompletion K _ z

theorem frobeniusCompletion_mem_int (K : PAdicLocalField p) (z : unramifiedCompletion K) :
    frobeniusCompletion K z ∈ unramifiedCompletionInt K ↔ z ∈ unramifiedCompletionInt K :=
  unramGalCompletion_mem_int K _ z

/-- **★★★★★★★★★★★★★★★★★★(Λ5)Frobenius は `K̂^{ur}` の連続な `K`-代数自己同型に延びる**。

`Gal(K^ur/K)` の位相的生成元 `σ`(`exists_coherentFrobenius`)は等長なので
一意に連続延長し、その延長は等長で `K` を固定する。

退化の自己検査。

* `Continuous Φ` を落とすと主張は弱くなりすぎる(`K^ur` 上の値だけでは
  `Φ` が決まらない)。実際 `unramGalCompletion_ext` は連続性を使っている。
* `∀ z, ‖Φ z‖ = ‖z‖` を落とすと `𝒪_{K̂^{ur}}` の保存が言えない。
* `σ` を `∃` の外に出すと**自由なパラメータ**になる(10 例目の退化と同型)ので、
  内側に閉じ込めてある。 -/
theorem exists_frobenius_lift (K : PAdicLocalField p) :
    ∃ Φ : unramifiedCompletion K ≃ₐ[K.carrier] unramifiedCompletion K,
      Continuous Φ ∧ (∀ z, ‖Φ z‖ = ‖z‖) ∧
      ∃ σ : unramGal K, (∀ N : ℕ, N ≠ 0 → σ ∈ unramLevelGeneratorSet K N) ∧
        ∀ x : ↥(unramifiedClosure K),
          Φ (x : unramifiedCompletion K)
            = ((σ x : ↥(unramifiedClosure K)) : unramifiedCompletion K) :=
  ⟨frobeniusCompletion K, continuous_frobeniusCompletion K, norm_frobeniusCompletion K,
    coherentFrobenius K, coherentFrobenius_spec K, frobeniusCompletion_coe K⟩

/-! ## 6. 剰余体——完備化は剰余体を変えない -/

/-- **`𝒪_{K^{ur}}`** —— 不分岐閉包の整数環。 -/
noncomputable def unramifiedClosureInt (K : PAdicLocalField p) :
    ValuationSubring ↥(unramifiedClosure K) :=
  (NormedField.valuation : Valuation ↥(unramifiedClosure K) NNReal).valuationSubring

@[simp] theorem mem_unramifiedClosureInt (K : PAdicLocalField p) (x : ↥(unramifiedClosure K)) :
    x ∈ unramifiedClosureInt K ↔ ‖x‖ ≤ 1 := by
  rw [unramifiedClosureInt, Valuation.mem_valuationSubring_iff, NormedField.valuation_apply]
  exact ⟨fun h => by exact_mod_cast h, fun h => by exact_mod_cast h⟩

theorem not_isUnit_unramifiedClosureInt (K : PAdicLocalField p) (w : ↥(unramifiedClosureInt K)) :
    ¬ IsUnit w ↔ ‖(w : ↥(unramifiedClosure K))‖ < 1 := by
  have h := Valuation.Integer.not_isUnit_iff_valuation_lt_one
    (v := (NormedField.valuation : Valuation ↥(unramifiedClosure K) NNReal)) (x := w)
  rw [NormedField.valuation_apply] at h
  exact h.trans ⟨fun hh => by exact_mod_cast hh, fun hh => by exact_mod_cast hh⟩

theorem not_isUnit_unramifiedCompletionInt (K : PAdicLocalField p)
    (w : ↥(unramifiedCompletionInt K)) :
    ¬ IsUnit w ↔ ‖(w : unramifiedCompletion K)‖ < 1 := by
  have h := Valuation.Integer.not_isUnit_iff_valuation_lt_one
    (v := (Valued.v : Valuation (unramifiedCompletion K) NNReal)) (x := w)
  rw [valuation_unramifiedCompletion] at h
  exact h.trans ⟨fun hh => by exact_mod_cast hh, fun hh => by exact_mod_cast hh⟩

theorem mem_maximalIdeal_unramifiedClosureInt (K : PAdicLocalField p)
    (w : ↥(unramifiedClosureInt K)) :
    w ∈ IsLocalRing.maximalIdeal ↥(unramifiedClosureInt K)
      ↔ ‖(w : ↥(unramifiedClosure K))‖ < 1 :=
  (IsLocalRing.mem_maximalIdeal w).trans (not_isUnit_unramifiedClosureInt K w)

theorem mem_maximalIdeal_unramifiedCompletionInt (K : PAdicLocalField p)
    (w : ↥(unramifiedCompletionInt K)) :
    w ∈ IsLocalRing.maximalIdeal ↥(unramifiedCompletionInt K)
      ↔ ‖(w : unramifiedCompletion K)‖ < 1 :=
  (IsLocalRing.mem_maximalIdeal w).trans (not_isUnit_unramifiedCompletionInt K w)

/-- `𝒪_{K^{ur}} → 𝒪_{K̂^{ur}}`。 -/
noncomputable def unramifiedIntHom (K : PAdicLocalField p) :
    ↥(unramifiedClosureInt K) →+* ↥(unramifiedCompletionInt K) :=
  RingHom.codRestrict
    ((UniformSpace.Completion.coeRingHom (α := ↥(unramifiedClosure K))).comp
      (SubringClass.subtype (unramifiedClosureInt K)))
    (unramifiedCompletionInt K)
    (fun x => by
      rw [mem_unramifiedCompletionInt]
      show ‖((x : ↥(unramifiedClosure K)) : unramifiedCompletion K)‖ ≤ 1
      rw [norm_coe_unramifiedCompletion]
      exact (mem_unramifiedClosureInt K _).mp x.2)

@[simp] theorem unramifiedIntHom_coe (K : PAdicLocalField p) (x : ↥(unramifiedClosureInt K)) :
    ((unramifiedIntHom K x : ↥(unramifiedCompletionInt K)) : unramifiedCompletion K)
      = ((x : ↥(unramifiedClosure K)) : unramifiedCompletion K) := rfl

instance isLocalHom_unramifiedIntHom (K : PAdicLocalField p) :
    IsLocalHom (unramifiedIntHom K) := by
  constructor
  intro a ha
  by_contra hcon
  rw [not_isUnit_unramifiedClosureInt] at hcon
  refine absurd ha ?_
  rw [not_isUnit_unramifiedCompletionInt, unramifiedIntHom_coe, norm_coe_unramifiedCompletion]
  exact hcon

/-- **剰余体の写像は全射**——`‖z‖ ≤ 1` なら `‖z - x‖ < 1` なる `x ∈ 𝒪_{K^{ur}}` が取れる。 -/
theorem residueField_map_surjective (K : PAdicLocalField p) :
    Function.Surjective (IsLocalRing.ResidueField.map (unramifiedIntHom K)) := by
  intro y
  obtain ⟨z, rfl⟩ := Ideal.Quotient.mk_surjective (I := IsLocalRing.maximalIdeal
    ↥(unramifiedCompletionInt K)) y
  obtain ⟨x, hx1, hx2⟩ := exists_norm_le_one_norm_sub_lt K (z : unramifiedCompletion K)
    ((mem_unramifiedCompletionInt K _).mp z.2) one_pos
  refine ⟨IsLocalRing.residue _ ⟨x, (mem_unramifiedClosureInt K x).mpr hx1⟩, ?_⟩
  rw [IsLocalRing.ResidueField.map_residue]
  refine Ideal.Quotient.eq.mpr ?_
  rw [mem_maximalIdeal_unramifiedCompletionInt]
  show ‖((unramifiedIntHom K ⟨x, _⟩ : ↥(unramifiedCompletionInt K)) : unramifiedCompletion K)
      - (z : unramifiedCompletion K)‖ < 1
  rw [unramifiedIntHom_coe, norm_sub_rev]
  exact hx2

/-- **★★★★★★★★★★★★★★(Λ5)`𝓀_{K^{ur}} ≅ 𝓀_{K̂^{ur}}`**——完備化は剰余体を変えない。

単射は「体からの環準同型は単射」、全射は `residueField_map_surjective`。

退化の自己検査:`K^ur` を有限次不分岐拡大 `K_N` に替えても主張は正しいが
(その場合はそもそも完備)、`K^ur` を `K.closure` に替えると
`𝓀_{K̂} = 𝔽̄_q` と `𝓀_{K.closure} = 𝔽̄_q` で依然正しい——本補題は
「代数拡大の完備化」に共通の現象であり、不分岐性は使っていない。 -/
noncomputable def residueFieldEquiv (K : PAdicLocalField p) :
    IsLocalRing.ResidueField ↥(unramifiedClosureInt K)
      ≃+* IsLocalRing.ResidueField ↥(unramifiedCompletionInt K) :=
  RingEquiv.ofBijective (IsLocalRing.ResidueField.map (unramifiedIntHom K))
    ⟨(IsLocalRing.ResidueField.map (unramifiedIntHom K)).injective,
      residueField_map_surjective K⟩

/-! ## 7. Frobenius の `𝒪_{K̂^{ur}}` への制限(Λ6 の舞台) -/

theorem norm_unramGalCompletion_symm (K : PAdicLocalField p) (σ : unramGal K)
    (z : unramifiedCompletion K) : ‖(unramGalCompletion K σ).symm z‖ = ‖z‖ := by
  conv_rhs => rw [← (unramGalCompletion K σ).apply_symm_apply z]
  rw [norm_unramGalCompletion]

/-- **Galois 作用の `𝒪_{K̂^{ur}}` への制限**。 -/
noncomputable def unramGalCompletionInt (K : PAdicLocalField p) (σ : unramGal K) :
    ↥(unramifiedCompletionInt K) ≃+* ↥(unramifiedCompletionInt K) where
  toFun w := ⟨unramGalCompletion K σ (w : unramifiedCompletion K),
    (unramGalCompletion_mem_int K σ _).mpr w.2⟩
  invFun w := ⟨(unramGalCompletion K σ).symm (w : unramifiedCompletion K), by
    rw [mem_unramifiedCompletionInt, norm_unramGalCompletion_symm]
    exact (mem_unramifiedCompletionInt K _).mp w.2⟩
  left_inv w := Subtype.ext ((unramGalCompletion K σ).symm_apply_apply _)
  right_inv w := Subtype.ext ((unramGalCompletion K σ).apply_symm_apply _)
  map_mul' a b := Subtype.ext (map_mul (unramGalCompletion K σ) _ _)
  map_add' a b := Subtype.ext (map_add (unramGalCompletion K σ) _ _)

@[simp] theorem unramGalCompletionInt_coe (K : PAdicLocalField p) (σ : unramGal K)
    (w : ↥(unramifiedCompletionInt K)) :
    ((unramGalCompletionInt K σ w : ↥(unramifiedCompletionInt K)) : unramifiedCompletion K)
      = unramGalCompletion K σ (w : unramifiedCompletion K) := rfl

/-- **`K̂^{ur}` の Frobenius の整数環への制限**。★Λ6(Dwork の補題)の舞台。

Dwork の補題は「`u ∈ 𝒪_{K̂^{ur}}^×` に対し `Φ(ξ)/ξ = u` なる `ξ ∈ 𝒪_{K̂^{ur}}^×` が
在る」という形で、本 `frobeniusCompletionInt` に関する主張になる。
逐次近似の材料(完備性・`exists_norm_le_one_norm_sub_lt`)は揃っているが、
各段の補正には剰余体上の Frobenius が `x ↦ x^q` であることが要る
(モジュール docstring「到達していない点」参照)。 -/
noncomputable def frobeniusCompletionInt (K : PAdicLocalField p) :
    ↥(unramifiedCompletionInt K) ≃+* ↥(unramifiedCompletionInt K) :=
  unramGalCompletionInt K (coherentFrobenius K)

end ABC3.Found.PGC
