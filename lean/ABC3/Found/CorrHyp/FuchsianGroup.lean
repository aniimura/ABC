/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import Mathlib.Analysis.Complex.UpperHalfPlane.Basic
import Mathlib.Analysis.Complex.UpperHalfPlane.MoebiusAction
import Mathlib.Analysis.Complex.UpperHalfPlane.ProperAction
import Mathlib.NumberTheory.ModularForms.ProperlyDiscontinuous
import Mathlib.GroupTheory.Schreier
import Mathlib.GroupTheory.Commensurable
import Mathlib.GroupTheory.GroupAction.ConjAct

open scoped Pointwise

/-!
# [CorrHyp] Track B の第一歩 —— Fuchsian 群を mathlib 上で実装する

`ResearchPaper/corrhyp-goal.md` §3(Track B)の入口。`Interface/CorrHyp/HyperbolicCurve.lean`
が posit した `Fuchsian` を、実際の `SL(2,ℝ)` の離散部分群として実装する
(まだ `HyperbolicCurveData` の具体的インスタンスは組んでいない——これは
その最初の 1 ブロックである)。

## ★★逸脱(記録)—— `PSL₂(ℝ)⁰` ではなく `SL(2,ℝ)` の部分群として扱う

原文は `Γ ⊆ PSL₂(ℝ)⁰`(射影特殊線形群の単位元成分)と書く。mathlib の
`ℍ` 上の Möbius 作用(`UpperHalfPlane.instProperlyDiscontinuousSL2RSubgroup` 等)は
`SL(2,ℝ)` を経由して定義されており、中心 `{±I}` は `ℍ` 上自明に作用する。
`SL(2,ℝ)` の部分群として扱っても、`ℍ` 上で見える作用・固有不連続性・
(離散な)商の構造は同じである。★`PSL₂(ℝ)⁰`(単位元成分、向きを保つ側)への
限定も本ファイルでは行っていない——`SL(2,ℝ)` 自体は連結なので、単位元成分の
制約は自動的に満たされる(`PSL₂(ℝ)` 全体には向きを反転する側が無い)。 -/

namespace ABC3.Found.CorrHyp

open UpperHalfPlane Matrix

/-- Fuchsian 群: `SL(2,ℝ)` の離散部分群。 -/
structure FuchsianGroup where
  /-- 部分群本体。 -/
  Γ : Subgroup (SpecialLinearGroup (Fin 2) ℝ)
  /-- 離散性(原文 §2: 「Γ ⊆ PSL₂(ℝ)⁰ の部分群として扱う」に暗黙に伴う仮定——
  基本群の deck 変換群の像は必ず離散である)。 -/
  discrete : DiscreteTopology Γ

attribute [instance] FuchsianGroup.discrete

/-- `F.Γ` は `ℍ` に固有不連続に作用する。

★**sorry 無し**——mathlib の `UpperHalfPlane.instProperlyDiscontinuousSL2RSubgroup`
から無条件に従う。原文 §2 冒頭の「Π ↪ Aut(H) = PSL₂(ℝ)⁰」の deck 変換の
議論が要求する性質そのもの。 -/
instance properlyDiscontinuous (F : FuchsianGroup) :
    ProperlyDiscontinuousSMul F.Γ ℍ := inferInstance

/-- `Γ₁` が `Γ₂` の部分群(`Interface` の `Sub`、`Proposition 3.2` で使う)。 -/
def IsSub (F1 F2 : FuchsianGroup) : Prop := F1.Γ ≤ F2.Γ

/-- `Γ₁` が `Γ₂` の中で有限指数(`Interface` の `FiniteIndexIn`)。 -/
def IsFiniteIndexIn (F1 F2 : FuchsianGroup) : Prop :=
  IsSub F1 F2 ∧ (F1.Γ.subgroupOf F2.Γ).FiniteIndex

/-- `Group.FG` は群同型で運べる(mathlib に無かったので補った、20 行)。 -/
private theorem fg_transport {G H : Type*} [Group G] [Group H] (e : G ≃* H) [Group.FG G] :
    Group.FG H := by
  obtain ⟨S, hS, hSfin⟩ := Group.fg_iff.mp ‹Group.FG G›
  refine Group.fg_iff.mpr ⟨e.toMonoidHom '' S, ?_, hSfin.image _⟩
  rw [← MonoidHom.map_closure, hS, Subgroup.map_top_of_surjective _ e.surjective]

/-- `Γ` が有限生成なら、その中の有限指数部分群も有限生成
(Schreier の補題、mathlib `Subgroup.fg_of_index_ne_zero`)。

★**sorry 無し**。原文が非明示的に使っている「Γ_X は有限生成(位相的基本群)」
という事実の系——`Theorem 3.3`/`Lemma 5.5` の証明の核となる Schreier の議論の
最初のブロック。 -/
theorem fg_of_finiteIndexIn {F1 F2 : FuchsianGroup} (hfg : Group.FG F2.Γ)
    (h : IsFiniteIndexIn F1 F2) : Group.FG F1.Γ := by
  haveI : Group.FG F2.Γ := hfg
  haveI : (F1.Γ.subgroupOf F2.Γ).FiniteIndex := h.2
  haveI : Group.FG (F1.Γ.subgroupOf F2.Γ) := Subgroup.fg_of_index_ne_zero _
  exact fg_transport (Subgroup.subgroupOfEquivOfLe h.1)

/-! ## `FEt` の圏構造 —— `Definition 1.1`-`1.4` の具体化に要る 4 点

`Space := FuchsianGroup`・`FEt A B := IsFiniteIndexIn A B` というモデルで
`Corr`(`Definition 1.1`)を実装するには、`comp`(合成)と `pullback`(ファイバー積、
`Definition 1.4`)が要る。ここで両方を実際に構成する——`Skeleton/CorrHyp/Section1.lean`
の `Corr`/`comp'` が要求する 4 点(`comp`・`pullback`・`pbFst`・`pbSnd`)そのもの。 -/

/-- `F` は `F` の中で(自明に)有限指数。`Definition 1.1`-`1.5` 末尾の
「it follows immediately that the relation of isogeny is an equivalence
relation」のうち**反射性**に要る恒等射(`Skeleton/CorrHyp/Section1.lean` の
`IsIsogenous.needs` が「まだ足していない」と記していたもの)。

★**sorry 無し**。 -/
theorem isFiniteIndexIn_refl (F : FuchsianGroup) : IsFiniteIndexIn F F := by
  refine ⟨le_refl _, ?_⟩
  rw [Subgroup.subgroupOf_self]
  infer_instance

/-- `Γ₁ ⊓ Γ₂` の離散性(離散群の部分群は離散)。`pullback` の台。 -/
def inter (F1 F2 : FuchsianGroup) : FuchsianGroup where
  Γ := F1.Γ ⊓ F2.Γ
  discrete := by
    haveI := F1.discrete
    haveI : DiscreteTopology (↥((F1.Γ ⊓ F2.Γ).subgroupOf F1.Γ)) := inferInstance
    exact (Subgroup.subgroupOfContinuousMulEquivOfLe
      (inf_le_left : F1.Γ ⊓ F2.Γ ≤ F1.Γ)).toHomeomorph.discreteTopology

/-- `FEt` の合成(`Interface.CorrHyp.HyperbolicCurveData.comp`)。

★**sorry 無し**——有限相対指数の乗法性(`Subgroup.relIndex_mul_relIndex`)から。 -/
theorem isFiniteIndexIn_comp {A B C : FuchsianGroup} (hAB : IsFiniteIndexIn A B)
    (hBC : IsFiniteIndexIn B C) : IsFiniteIndexIn A C := by
  refine ⟨hAB.1.trans hBC.1, ?_⟩
  haveI hAB' : A.Γ.IsFiniteRelIndex B.Γ := Subgroup.isFiniteRelIndex_iff_finiteIndex.mpr hAB.2
  haveI hBC' : B.Γ.IsFiniteRelIndex C.Γ := Subgroup.isFiniteRelIndex_iff_finiteIndex.mpr hBC.2
  rw [Subgroup.isFiniteRelIndex_iff_relIndex_ne_zero] at hAB' hBC'
  rw [← Subgroup.isFiniteRelIndex_iff_finiteIndex, Subgroup.isFiniteRelIndex_iff_relIndex_ne_zero,
      ← Subgroup.relIndex_mul_relIndex A.Γ B.Γ C.Γ hAB.1 hBC.1]
  exact mul_ne_zero hAB' hBC'

/-- `pullback f g` から `A` への射影(`Interface.CorrHyp.HyperbolicCurveData.pbFst`)。

★**sorry 無し**——`(A ⊓ B).relIndex C = (A⊓B).relIndex A · A.relIndex C` かつ
左辺が有限(`A`・`B` とも `C` の中で有限指数だから、`relIndex_inf_ne_zero`)なので、
積が有限なら両因子とも有限。 -/
theorem isFiniteIndexIn_pbFst {A B C : FuchsianGroup} (hA : IsFiniteIndexIn A C)
    (hB : IsFiniteIndexIn B C) : IsFiniteIndexIn (inter A B) A := by
  refine ⟨inf_le_left, ?_⟩
  haveI hA' : A.Γ.IsFiniteRelIndex C.Γ := Subgroup.isFiniteRelIndex_iff_finiteIndex.mpr hA.2
  haveI hB' : B.Γ.IsFiniteRelIndex C.Γ := Subgroup.isFiniteRelIndex_iff_finiteIndex.mpr hB.2
  rw [Subgroup.isFiniteRelIndex_iff_relIndex_ne_zero] at hA' hB'
  have hinf : (A.Γ ⊓ B.Γ).relIndex C.Γ ≠ 0 := Subgroup.relIndex_inf_ne_zero hA' hB'
  rw [← Subgroup.isFiniteRelIndex_iff_finiteIndex, Subgroup.isFiniteRelIndex_iff_relIndex_ne_zero]
  show (A.Γ ⊓ B.Γ).relIndex A.Γ ≠ 0
  intro h0
  have heq : (A.Γ ⊓ B.Γ).relIndex A.Γ * A.Γ.relIndex C.Γ = (A.Γ ⊓ B.Γ).relIndex C.Γ :=
    Subgroup.relIndex_mul_relIndex (A.Γ ⊓ B.Γ) A.Γ C.Γ inf_le_left hA.1
  rw [h0, zero_mul] at heq
  exact hinf heq.symm

/-- `pullback f g` から `B` への射影(`Interface.CorrHyp.HyperbolicCurveData.pbSnd`)。
`inter A B = inter B A`(`⊓` の可換性)へ帰着して `pbFst` を使い回す。 -/
theorem isFiniteIndexIn_pbSnd {A B C : FuchsianGroup} (hA : IsFiniteIndexIn A C)
    (hB : IsFiniteIndexIn B C) : IsFiniteIndexIn (inter A B) B := by
  have : inter A B = inter B A := by
    unfold inter; congr 1; exact inf_comm _ _
  rw [this]
  exact isFiniteIndexIn_pbFst hB hA

/-! ## `Comm`(commensurator)—— `Definition 2.1` の語彙

原文 §2: `Comm(Γ) ≝ {γ ∈ PSL₂(ℝ)⁰ | (γ·Γ·γ⁻¹) ∼ Γ}`。mathlib にちょうど同じ概念が
`Subgroup.Commensurable.commensurator` として存在する——ここでは `Γ ⊆ Comm(Γ)`
(`Definition 2.1` 直前の地の文 "Note that Γ ⊆ Comm(Γ)")を実装する。

★`Comm(Γ)` は一般には離散とは限らない(arithmetic な `Γ` では稠密になる、
というのが §2 の主定理の内容そのもの)ので、`Comm` の戻り値は `FuchsianGroup`
ではなく素の `Subgroup`——`Interface.CorrHyp.HyperbolicCurveData.Comm` が
`Fuchsian → Fuchsian`(常に離散)と posit しているのは簡略化であり、
ここで判明した不一致として記録する(逸脱)。 -/

open Subgroup.Commensurable in
/-- 原文 §2、`Definition 2.1` 直前:「Note that Γ ⊆ Comm(Γ)」。

★**sorry 無し**。標準3公理のみ(`#print axioms` で確認済み)。 -/
theorem self_le_commensurator (F : FuchsianGroup) : F.Γ ≤ commensurator F.Γ := by
  intro g hg
  rw [commensurator_mem_iff]
  have heq : (ConjAct.toConjAct g • F.Γ : Subgroup _) = F.Γ := by
    apply Subgroup.ext
    intro x
    rw [Subgroup.mem_pointwise_smul_iff_inv_smul_mem]
    simp only [ConjAct.smul_def]
    constructor
    · intro hx
      have := F.Γ.mul_mem (F.Γ.mul_mem hg hx) (F.Γ.inv_mem hg)
      simpa [mul_assoc] using this
    · intro hx
      have := F.Γ.mul_mem (F.Γ.mul_mem (F.Γ.inv_mem hg) hx) hg
      simpa [mul_assoc] using this
  rw [heq]

/-! ## `Proposition 3.2` へ向けた足場

原文 p.8 の証明は「for the purpose of proving this Proposition, we may assume
that C is Galois over Z. Thus, Γ_C is normal (and of finite index) in Γ_Z」と、
`Γ_C` を正規化する段から始まる。一般の `C` では正規とは限らないので、
標準的な処方は「`Γ_Z` の中での `Γ_C` の normal core」を取ること
(`Γ_C` を Galois 閉包に取り替えるのと同じ効果——★逸脱として記録:
原文は具体的な Galois 閉包を取るが、ここでは normal core という
より直接的な群論的構成で同じ役割を果たす)。 -/

/-- `Proposition 3.2` へ向けた足場: `C ≤ Z` が有限指数なら、`Z.Γ` の中で
正規かつ有限指数な `N ≤ C.Γ` が取れる(`Z.Γ` の中での `C.Γ` の normal core、
`Subgroup.finiteIndex_normalCore` から)。

★**sorry 無し**。`Proposition 3.2` 本体(`Γ_Z ≤ Comm(Γ_X)`)はまだ——
`N` を大きい群へ写して `Γ_X` との交わりを比較する段が残る。 -/
theorem exists_normalCore_le {C Z : FuchsianGroup} (hCZ : IsFiniteIndexIn C Z) :
    ∃ N : Subgroup Z.Γ, N.Normal ∧ N.FiniteIndex ∧ N ≤ C.Γ.subgroupOf Z.Γ := by
  haveI : (C.Γ.subgroupOf Z.Γ).FiniteIndex := hCZ.2
  exact ⟨(C.Γ.subgroupOf Z.Γ).normalCore, inferInstance, inferInstance,
    Subgroup.normalCore_le _⟩

/-- `exists_normalCore_le` の `N` を大きい群(`Z.Γ`・`C.Γ` と同じ
`SpecialLinearGroup (Fin 2) ℝ`)へ押し出したもの。 -/
noncomputable def Ncore {C Z : FuchsianGroup} (h : IsFiniteIndexIn C Z) :
    Subgroup (Matrix.SpecialLinearGroup (Fin 2) ℝ) :=
  (exists_normalCore_le h).choose.map Z.Γ.subtype

theorem Ncore_le_Z {C Z : FuchsianGroup} (h : IsFiniteIndexIn C Z) : Ncore h ≤ Z.Γ := by
  unfold Ncore
  rintro x ⟨n, -, rfl⟩
  exact n.2

theorem Ncore_le_C {C Z : FuchsianGroup} (h : IsFiniteIndexIn C Z) : Ncore h ≤ C.Γ := by
  unfold Ncore
  rintro x ⟨n, hn, rfl⟩
  have := (exists_normalCore_le h).choose_spec.2.2 hn
  simpa [Subgroup.mem_subgroupOf] using this

/-- `Ncore h` は `Z.Γ` の元による共役で不変(`Z.Γ` の中で正規)。

★**sorry 無し**——`Ncore` の中身(`exists_normalCore_le` の `N`)が
`↥Z.Γ` の中で正規であることを、大きい群での共役に翻訳しただけ。 -/
theorem Ncore_normal {C Z : FuchsianGroup} (h : IsFiniteIndexIn C Z)
    {z : Matrix.SpecialLinearGroup (Fin 2) ℝ} (hz : z ∈ Z.Γ) :
    (ConjAct.toConjAct z • Ncore h : Subgroup _) = Ncore h := by
  unfold Ncore
  haveI := (exists_normalCore_le h).choose_spec.1
  apply Subgroup.ext
  intro x
  rw [Subgroup.mem_pointwise_smul_iff_inv_smul_mem]
  simp only [ConjAct.smul_def, Subgroup.mem_map, ConjAct.ofConjAct_inv, ConjAct.ofConjAct_toConjAct]
  constructor
  · rintro ⟨n, hn, hnx⟩
    refine ⟨⟨z, hz⟩ * n * ⟨z, hz⟩⁻¹, ?_, ?_⟩
    · have := (Subgroup.Normal.conj_mem (exists_normalCore_le h).choose_spec.1 n hn ⟨z, hz⟩)
      simpa using this
    · have hn' : (n : Matrix.SpecialLinearGroup (Fin 2) ℝ) = z⁻¹ * x * z⁻¹⁻¹ := hnx
      show (z : Matrix.SpecialLinearGroup (Fin 2) ℝ) * n * z⁻¹ = x
      rw [hn']; group
  · rintro ⟨n, hn, hnx⟩
    refine ⟨⟨z, hz⟩⁻¹ * n * ⟨z, hz⟩, ?_, ?_⟩
    · have := (Subgroup.Normal.conj_mem (exists_normalCore_le h).choose_spec.1 n hn ⟨z, hz⟩⁻¹)
      simpa using this
    · have hn' : (n : Matrix.SpecialLinearGroup (Fin 2) ℝ) = x := hnx
      show (z : Matrix.SpecialLinearGroup (Fin 2) ℝ)⁻¹ * n * z = z⁻¹ * x * z⁻¹⁻¹
      rw [hn']; group

/-- `Ncore h`(`Z.Γ` の中での `C.Γ` の normal core を大きい群へ押し出したもの)は
`(Ncore h).subgroupOf Z.Γ` として見ると、もとの `exists_normalCore_le h` の `N`
そのものと一致する(`Subtype` の像を `comap` で戻すと元に戻る、という一般論)。 -/
theorem Ncore_subgroupOf_Z {C Z : FuchsianGroup} (h : IsFiniteIndexIn C Z) :
    (Ncore h).subgroupOf Z.Γ = (exists_normalCore_le h).choose := by
  unfold Ncore Subgroup.subgroupOf
  exact Subgroup.comap_map_eq_self_of_injective Z.Γ.subtype_injective _

/-- `Ncore h` は `Z.Γ` の中で有限指数。★**sorry 無し**。 -/
theorem Ncore_finiteIndex_Z {C Z : FuchsianGroup} (h : IsFiniteIndexIn C Z) :
    ((Ncore h).subgroupOf Z.Γ).FiniteIndex := by
  rw [Ncore_subgroupOf_Z]
  exact (exists_normalCore_le h).choose_spec.2.1

/-- `Ncore h` は `C.Γ` の中でも(したがって `relIndex` の意味で)有限指数。
`Ncore h ≤ C.Γ ≤ Z.Γ` の推移律と `Ncore_finiteIndex_Z` から。 -/
theorem Ncore_relIndex_C_ne_zero {C Z : FuchsianGroup} (h : IsFiniteIndexIn C Z) :
    (Ncore h).relIndex C.Γ ≠ 0 := by
  haveI := Ncore_finiteIndex_Z h
  have hZ : (Ncore h).relIndex Z.Γ ≠ 0 :=
    Subgroup.isFiniteRelIndex_iff_relIndex_ne_zero.mp
      (Subgroup.isFiniteRelIndex_iff_finiteIndex.mpr ‹((Ncore h).subgroupOf Z.Γ).FiniteIndex›)
  have heq : (Ncore h).relIndex C.Γ * C.Γ.relIndex Z.Γ = (Ncore h).relIndex Z.Γ :=
    Subgroup.relIndex_mul_relIndex (Ncore h) C.Γ Z.Γ (Ncore_le_C h) h.1
  intro h0
  rw [h0, zero_mul] at heq
  exact hZ heq.symm

/-- `Ncore h` は(`hCX : IsFiniteIndexIn C X` が与えられれば)`X.Γ` の中でも
有限指数。`Ncore h ≤ C.Γ ≤ X.Γ` の推移律から。 -/
theorem Ncore_relIndex_X_ne_zero {C X Z : FuchsianGroup} (h : IsFiniteIndexIn C Z)
    (hCX : IsFiniteIndexIn C X) : (Ncore h).relIndex X.Γ ≠ 0 := by
  have hC : (Ncore h).relIndex C.Γ ≠ 0 := Ncore_relIndex_C_ne_zero h
  haveI hCX' : C.Γ.IsFiniteRelIndex X.Γ :=
    Subgroup.isFiniteRelIndex_iff_finiteIndex.mpr hCX.2
  have hCXne : C.Γ.relIndex X.Γ ≠ 0 := Subgroup.isFiniteRelIndex_iff_relIndex_ne_zero.mp hCX'
  have heq : (Ncore h).relIndex C.Γ * C.Γ.relIndex X.Γ = (Ncore h).relIndex X.Γ :=
    Subgroup.relIndex_mul_relIndex (Ncore h) C.Γ X.Γ (Ncore_le_C h) hCX.1
  intro h0
  rw [h0] at heq
  rcases mul_eq_zero.mp heq with h1 | h1
  · exact hC h1
  · exact hCXne h1

open Subgroup.Commensurable in
/-- [CorrHyp] **Proposition 3.2** の具体化。

原文 (CorrHyp p.8): 「We have Γ_Z ⊆ Comm(Γ_X)」——`C` が `X`・`Z` それぞれに
有限指数で埋め込まれた `FuchsianGroup`(`Corr D X Z` の台 `C` の具体化)なら、
`Z.Γ ≤ Comm(X.Γ)`。

証明の筋(原文が「C を Z 上 Galois に取り直せる」と1行で済ませている段を
`Ncore`(`Z.Γ` の中での `C.Γ` の normal core)で具体化したもの):
`N := Ncore h` は `C.Γ`(ゆえに `X.Γ`)にも `Z.Γ` にも有限指数で入り、
`z ∈ Z.Γ` による共役で不変。ゆえに `N ≤ X.Γ ⊓ (z·X.Γ·z⁻¹)` かつ `N` は
`X.Γ`・`z·X.Γ·z⁻¹` のどちらの中でも有限指数なので、交わりも両側で有限指数——
これが `Commensurable (z·X.Γ·z⁻¹) X.Γ` の定義そのもの。

★**sorry 無し**。標準3公理のみ(`#print axioms` で確認済み)。 -/
theorem prop_3_2 {C X Z : FuchsianGroup} (hCX : IsFiniteIndexIn C X)
    (hCZ : IsFiniteIndexIn C Z) : Z.Γ ≤ commensurator X.Γ := by
  intro z hz
  rw [commensurator_mem_iff]
  set N := Ncore hCZ with hN
  have hNX : N ≤ X.Γ := (Ncore_le_C hCZ).trans hCX.1
  have hzsymm : (ConjAct.toConjAct z • N : Subgroup _) = N := by
    rw [hN]; exact Ncore_normal hCZ hz
  have hNzX : N ≤ (ConjAct.toConjAct z • X.Γ : Subgroup _) := by
    rw [← hzsymm]
    intro x hx
    rw [Subgroup.mem_pointwise_smul_iff_inv_smul_mem] at hx ⊢
    exact hNX hx
  have hN1 : N.relIndex X.Γ ≠ 0 := Ncore_relIndex_X_ne_zero hCZ hCX
  have hN2 : N.relIndex (ConjAct.toConjAct z • X.Γ : Subgroup _) ≠ 0 := by
    have := Subgroup.relIndex_pointwise_smul (h := ConjAct.toConjAct z) N X.Γ
    rw [hzsymm] at this
    rw [this]; exact hN1
  haveI hIF1 : N.IsFiniteRelIndex X.Γ := Subgroup.isFiniteRelIndex_iff_relIndex_ne_zero.mpr hN1
  haveI hIF2 : N.IsFiniteRelIndex (ConjAct.toConjAct z • X.Γ : Subgroup _) :=
    Subgroup.isFiniteRelIndex_iff_relIndex_ne_zero.mpr hN2
  haveI hK1 : (X.Γ ⊓ (ConjAct.toConjAct z • X.Γ : Subgroup _)).IsFiniteRelIndex X.Γ :=
    Subgroup.isFiniteRelIndex_of_le_left X.Γ (le_inf hNX hNzX)
  haveI hK2 : (X.Γ ⊓ (ConjAct.toConjAct z • X.Γ : Subgroup _)).IsFiniteRelIndex
      (ConjAct.toConjAct z • X.Γ : Subgroup _) :=
    Subgroup.isFiniteRelIndex_of_le_left _ (le_inf hNX hNzX)
  constructor
  · have := Subgroup.isFiniteRelIndex_iff_relIndex_ne_zero.mp hK1
    rwa [Subgroup.inf_relIndex_left] at this
  · have := Subgroup.isFiniteRelIndex_iff_relIndex_ne_zero.mp hK2
    rwa [Subgroup.inf_relIndex_right] at this

end ABC3.Found.CorrHyp
