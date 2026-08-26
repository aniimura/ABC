/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Def24ScTransport
import ABC3.Found.FrdI.Prop55Rlf
import ABC3.Found.FrdI.Prop55Std

/-!
# [FrdI] Proposition 5.5, (iii) —— `𝒞^rlf` は rationally standard 型

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.105。

原文 (FrdI p.105):
> if C is of standard (respectively, rationally standard) type, then so are Cun-tr, Crlf.

## ★★★`Theorem 6.4, (i)` の 5 圏の**最後の 1 つ**

`Theorem 6.4, (i)` は `𝒞, 𝒞^pf, 𝒞^rlf, 𝒞^un-tr, (𝒞^pf)^un-tr` が
rationally standard だと言う。★他の 4 つは

| 圏 | 一般形 |
|---|---|
| `𝒞` | 模型ごと(`ModelData.model_isOfRationallyStandardType_of_baseId`) |
| `𝒞^un-tr` | `unTr_isOfRationallyStandardType_gen`(`Prop55UnTrRatStd.lean`) |
| `𝒞^pf` | `pfRoot_isOfRationallyStandardType`(`Prop55PfRatStd.lean`) |
| `(𝒞^pf)^un-tr` | `pfRoot_unTr_isOfRationallyStandardType`(`Prop55UnTrRatStd.lean`) |

で済んでおり、**本ファイルが `𝒞^rlf` を埋める**。

## ★★4 条の出どころ

`𝒞^rlf` は model Frobenioid **そのもの**(`crlf_eq_model`)なので、
`Definition 4.5, (iii)` の 4 条は次のように埋まる:

| 条 | 出どころ |
|---|---|
| (a) birat-Frobenius-normalized | `ModelData.model_isOfBiratFrobNormalizedType`(在庫) |
| (a) rational | ★`scModel_isOfStrictlyRationalType`(本ファイル) |
| (a) standard | `scModel_standardType`(在庫) |
| (b) Frobenius-compact 対象 | 仮引数(模型ごと) |

## ★★★rational が通った理由 —— 素点と台が `toSc` で移るから

`IsStrictlyRational` は**すべての素点**について量化するので、
`Prime(Φ^rlf)` と `Prime(Φ)` の対応が要る。★`𝒞^pf` のときの
`primeToPf_bijective` に当たるものが実化の側には無かったが、
`Def24ScTransport.lean` で

    primeToSc_bijective       素点は `toSc` で 1 対 1
    mem_suppElt_primeToSc     台も `toSc` でちょうど移る

を作ったので、`𝒞` の側の `hsp` をそのまま押し出せる。
★`Φ^birat` の側は `mem_sPhiBiratOn_of_phiBiratOn`(在庫)がそのまま効く
(`ℝ·Φ^birat` は `Φ^birat` の像を含むので)。

★★`ι^rlf` は**指示関数**(`iotaInd`)に取ってよい ——
`SuppElt` は `ι` を含まない条件と同値(`mem_supp_factorMap_iff`)だからである。
-/

namespace ABC3.Found.FrdI

open CategoryTheory
open scoped NNReal TensorProduct

universe v u w u2 v2

section RlfRatStd

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} {P : PreFrobenioid C Φ}

set_option maxHeartbeats 1000000 in
/-- ★★★★★★**`𝒞^rlf` は strictly rational 型**。

★`𝒞` の側の `hsp`(素点ごとの「正・負の分割」)を `toSc` で押し出すだけである:

| 移すもの | 使う補題 |
|---|---|
| 素点 | `primeToSc_bijective`(全射性) |
| 台 | `mem_suppElt_primeToSc` |
| `Φ^birat` | `mem_sPhiBiratOn_of_phiBiratOn`(在庫) |
-/
theorem scModel_isOfStrictlyRationalType (G : Frobenioid P)
    (hiso : ∀ Y : C, IsIsotropic P Y)
    (hfn : ∀ X : BiratCat P G, IsFrobeniusNormalized (biratPre P G) X)
    (hcharInj : ∀ {A B : D} (α : B ⟶ A),
      IsCharacteristicallyInjective (scMap (S := ℝ≥0) (Φ.map α)))
    (hint : ∀ A : D, IsIntegralMonoid (ScT ℝ≥0 (Φ.val A)))
    (hfsmD : IsOfFSMType D)
    (hdivRlf : (phiScOn ℝ≥0 Φ hcharInj).IsDivisorialOn)
    (htot : IsTotallyEpimorphic D) (hconn : IsConnected D)
    (ι : ∀ A : D, Prime (Φ.val A) → Pf (Φ.val A) → ℝ≥0)
    (Hpf : ∀ A : D, IsPerfFactorialWith (Φ.val A) (ι A))
    (hperf : ∀ A : D, IsPerfectMonoid (Φ.val A))
    (hdiv : ∀ A : D, IsDivisorial (Φ.val A))
    (hsp : ∀ (d : D) (p : Prime (Φ.val d)),
      ∃ (a b : Φ.val d) (y : Gp (Φ.val d)), y ∈ phiBiratOn G d ∧
        (toGp _ a - toGp _ b = y ∨ toGp _ a - toGp _ b = -y) ∧
        p ∈ SuppElt (ι d) a ∧ p ∉ SuppElt (ι d) b) :
    IsOfStrictlyRationalType
      (ModelData.Obj (scModel ℝ≥0 G hiso hfn hcharInj hint hfsmD))
      (ModelData.modelPre (scModel_hyp G hiso hfn hcharInj hint hfsmD hdivRlf htot hconn))
      (ModelData.model_frobenioid
        (scModel_hyp G hiso hfn hcharInj hint hfsmD hdivRlf htot hconn))
      (fun d => iotaInd (ScT ℝ≥0 (Φ.val d))) := by
  haveI := hconn
  set hyp := scModel_hyp G hiso hfn hcharInj hint hfsmD hdivRlf htot hconn with hhyp
  intro X q
  obtain ⟨p, hp⟩ := (primeToSc_bijective (Hpf X.base) (hperf X.base) (hdiv X.base)).2 q
  obtain ⟨a, b, y, hy, hab, hpa, hpb⟩ := hsp X.base p
  have hkey : ∀ z : Gp (Φ.val X.base), z ∈ phiBiratOn G X.base →
      toScGp (S := ℝ≥0) z ∈ phiBiratAt (ModelData.modelPre hyp)
        (ModelData.model_frobenioid hyp) (show BiratCat _ _ from X) := by
    intro z hz
    exact divB_mem_phiBiratAt (ModelData.modelRatFnData hyp)
      (A := X) ⟨toScGp z, mem_sPhiBiratOn_of_phiBiratOn G hz⟩
  have hsub : toGp (ScT ℝ≥0 (Φ.val X.base)) (toSc a) - toGp _ (toSc b)
      = toScGp (S := ℝ≥0) (toGp (Φ.val X.base) a - toGp _ b) := by
    show _ = gpMap _ toSc (toGp (Φ.val X.base) a - toGp _ b)
    rw [map_sub, gpMap_toGp, gpMap_toGp]
  have hmem : toScGp (S := ℝ≥0) (toGp (Φ.val X.base) a - toGp _ b)
      ∈ phiBiratAt (ModelData.modelPre hyp)
        (ModelData.model_frobenioid hyp) (show BiratCat _ _ from X) := by
    rcases hab with h | h
    · rw [h]; exact hkey y hy
    · rw [h, map_neg]
      exact AddSubgroup.neg_mem _ (hkey y hy)
  refine ⟨toSc a, toSc b, ?_, ?_, ?_⟩
  · exact Eq.mpr (congrArg (fun t => t ∈ phiBiratAt (ModelData.modelPre hyp)
      (ModelData.model_frobenioid hyp) (show BiratCat _ _ from X)) hsub) hmem
  · rw [← hp]
    exact (mem_suppElt_primeToSc (Hpf X.base) (hperf X.base) (hdiv X.base) a p).mpr hpa
  · rw [← hp]
    exact fun hc =>
      hpb ((mem_suppElt_primeToSc (Hpf X.base) (hperf X.base) (hdiv X.base) b p).mp hc)

set_option maxHeartbeats 1000000 in
/-- ★★★★★★★**[FrdI] Proposition 5.5, (iii)** —— **`𝒞^rlf` は rationally standard 型**。

★★これで `Theorem 6.4, (i)` が並べる 5 圏の**一般形がすべて揃った**
(`𝒞` は模型ごと、`𝒞^un-tr` / `𝒞^pf` / `(𝒞^pf)^un-tr` / `𝒞^rlf` は一般形)。

★(b)(`((𝒞^rlf)^un-tr)^birat` の Frobenius-compact 対象)だけは仮引数で受ける ——
在庫の `unTr_isOfRationallyStandardType` と同じ流儀で、模型ごとに与えるものである。

## ★★★★訂正(2026-08-25)—— `Skeletal 𝒟` は要らない

当初 (a) を `scModel_isOfModelType` の第 2 成分から取っていたが、
それは **`Skeletal 𝒟` を要求する**。★算術の `𝒟 = (FinSub F K̄)ᵒᵖ` では
**共役な中間体が同型だが相異なる**ので `Skeletal 𝒟` は**偽**であり、
それを仮定に置くと定理が空虚になる。

★`biratFrobNormalized` は `ModelData.model_isOfBiratFrobNormalizedType` から
**`Skeletal` 無しで**取れるので、そちらに差し替えた。 -/
theorem scModel_isOfRationallyStandardType (G : Frobenioid P)
    (hiso : ∀ Y : C, IsIsotropic P Y)
    (hfn : ∀ X : BiratCat P G, IsFrobeniusNormalized (biratPre P G) X)
    (hcharInj : ∀ {A B : D} (α : B ⟶ A),
      IsCharacteristicallyInjective (scMap (S := ℝ≥0) (Φ.map α)))
    (hint : ∀ A : D, IsIntegralMonoid (ScT ℝ≥0 (Φ.val A)))
    (hfsmD : IsOfFSMType D)
    (hdivRlf : (phiScOn ℝ≥0 Φ hcharInj).IsDivisorialOn)
    (htot : IsTotallyEpimorphic D) (hconn : IsConnected D)
    (hfsmff : IsOfFSMFFType D)
    (hnd : (scModel ℝ≥0 G hiso hfn hcharInj hint hfsmD).phi.IsNonDilatingOn)
    (hngl : ¬ IsOfGroupLikeType
      (ModelData.modelPre (scModel_hyp G hiso hfn hcharInj hint hfsmD hdivRlf htot hconn)))
    (ι : ∀ A : D, Prime (Φ.val A) → Pf (Φ.val A) → ℝ≥0)
    (Hpf : ∀ A : D, IsPerfFactorialWith (Φ.val A) (ι A))
    (hperf : ∀ A : D, IsPerfectMonoid (Φ.val A))
    (hdiv : ∀ A : D, IsDivisorial (Φ.val A))
    (hsp : ∀ (d : D) (p : Prime (Φ.val d)),
      ∃ (a b : Φ.val d) (y : Gp (Φ.val d)), y ∈ phiBiratOn G d ∧
        (toGp _ a - toGp _ b = y ∨ toGp _ a - toGp _ b = -y) ∧
        p ∈ SuppElt (ι d) a ∧ p ∉ SuppElt (ι d) b)
    (hcpt : ∃ X : BiratCat
        (unTrPre (ModelData.modelPre
            (scModel_hyp G hiso hfn hcharInj hint hfsmD hdivRlf htot hconn))
          (ModelData.model_frobenioid
            (scModel_hyp G hiso hfn hcharInj hint hfsmD hdivRlf htot hconn)).core)
        (unTr_frobenioid _ _ (ModelData.model_frobenioid
          (scModel_hyp G hiso hfn hcharInj hint hfsmD hdivRlf htot hconn))),
      IsFrobeniusCompact (unTrBiratPre _ _ (ModelData.model_frobenioid
        (scModel_hyp G hiso hfn hcharInj hint hfsmD hdivRlf htot hconn))) X) :
    haveI := hconn
    IsOfRationallyStandardType
      (ModelData.modelPre (scModel_hyp G hiso hfn hcharInj hint hfsmD hdivRlf htot hconn))
      (ModelData.model_frobenioid
        (scModel_hyp G hiso hfn hcharInj hint hfsmD hdivRlf htot hconn))
      (fun d => iotaInd (ScT ℝ≥0 (Φ.val d))) := by
  haveI := hconn
  exact
    { biratFrobNormalized :=
        ModelData.model_isOfBiratFrobNormalizedType
          (scModel_hyp G hiso hfn hcharInj hint hfsmD hdivRlf htot hconn)
      rational := fun X => isRational_of_isStrictlyRational _ X
        (scModel_isOfStrictlyRationalType G hiso hfn hcharInj hint hfsmD hdivRlf htot hconn
          ι Hpf hperf hdiv hsp X)
      standard := scModel_standardType G hiso hfn hcharInj hint hfsmD hdivRlf htot hconn
        _ hfsmff hnd hngl
      unTrBiratCompact := hcpt }

end RlfRatStd

/-! ### ★出典の紐付け -/

def scModel_isOfStrictlyRationalType.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 105,
    item := "Proposition 5.5, (iii) — 𝒞^rlf は strictly rational 型",
    sectionId := "frdi-prop-5-5" }

def scModel_isOfRationallyStandardType.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 105,
    item := "Proposition 5.5, (iii) — 𝒞^rlf は rationally standard 型",
    sectionId := "frdi-prop-5-5" }

def scModel_isOfRationallyStandardType.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "primeToSc_bijective(素点は toSc で 1 対 1)"
      (.inProject "ABC3" "ABC3.Found.FrdI.primeToSc_bijective") 12,
    .citation "[ABC3]" "mem_suppElt_primeToSc(台も toSc でちょうど移る)"
      (.inProject "ABC3" "ABC3.Found.FrdI.mem_suppElt_primeToSc") 48,
    .citation "[ABC3]" "scModel_isOfModelType((a) birat-Frobenius-normalized)"
      (.inProject "ABC3" "ABC3.Found.FrdI.scModel_isOfModelType") 104,
    .citation "[ABC3]" "scModel_standardType((a) standard)"
      (.inProject "ABC3" "ABC3.Found.FrdI.scModel_standardType") 105,
    .citation "[ABC3]" "mem_sPhiBiratOn_of_phiBiratOn(Φ^birat の像は ℝ·Φ^birat に入る)"
      (.inProject "ABC3" "ABC3.Found.FrdI.mem_sPhiBiratOn_of_phiBiratOn") 103,
    .implicitStep
      ("★(b)(Frobenius-compact 対象)は仮引数で受ける。" ++
       "また scModel の仮定束(hcharInj / Φ^rlf が divisorial / Skeletal 𝒟)も" ++
       "仮引数であり、模型ごとに与える(節点 rlf-flat と同じ宿題)") 105 ]

end ABC3.Found.FrdI
