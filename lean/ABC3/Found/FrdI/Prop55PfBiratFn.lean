/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Prop55ScaleRootBirat
import ABC3.Found.FrdI.Prop55Std
import ABC3.Found.FrdI.Prop55PfSlice

/-!
# [FrdI] Proposition 5.5, (iii) —— `𝒞^pf` は birationally Frobenius-normalized 型

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.105。

原文 (FrdI p.105):
> if C is of standard (respectively, rationally standard) type, then so are Cun-tr, Crlf.

## ★★このファイルの位置づけ —— `Theorem 6.4, (i)` の律速 1 点だった条

`Theorem 6.4, (i)` は `C, C^pf, C^rlf, C^un-tr, (C^pf)^un-tr` の 5 圏が
rationally standard だと言う。★`𝒞^pf` の `Definition 4.5, (iii)` の 4 条のうち
**standard**(`pfRoot_standardType`)・
**`((𝒞^pf)^un-tr)^birat` の Frobenius-compact 対象**(`unTr_pf_biratCompact_of_baseId`)・
**rational の台の対応**(`Def24PfTransport.lean`)は先に済んでおり、
最後に残ったのが**本ファイルの birationally Frobenius-normalized** である。

## ★★★なぜ在庫の `pfRoot_isOfModelType` では出ないか

在庫の `pfRoot_isOfModelType_of_arb` は 4 条を要求し、そのうち
**unit-trivial 型**だけが `Example 6.1` / `Example 6.3` の `𝒞` で成り立たない
(単元が `k_L^×` / `μ(L)`)。★これが `hut` の壁であり、本ファイルはそこを迂回する。

## ★★★★証明の骨格(3 段 ＋ 移送)

★**`isFrobeniusNormalized_transport`(在庫)は関手を要求しない** ——
1 つの対象での `End` の乗法同型 ＋ 3 つの対応だけである。

| 段 | 写像 | 乗法性の根拠 | `Base`・`degFr` の保存 |
|---|---|---|---|
| 1 | `biratPfHom` | `biratPfHom_comp` | `biratBase_biratPfHom_eq` / `biratDeg_biratPfHom_eq` |
| 2 | `scaleRootBirat (n·n)` | **関手**(`Functor.mapEnd`) | `biratBase_scaleRootBirat` / `biratDeg_scaleRootBirat` |
| 3 | 同型 `⟨A,n⟩ ≅ ⟨A^{(n)}, n·n⟩` による共役 | —— | `isFrobeniusNormalized_of_iso` |

★1・2 段を `pfBiratEndMulEquiv12` に束ねて `isFrobeniusNormalized_transport` に流し、
3 段目は在庫の `isFrobeniusNormalized_of_iso` で済ませる
(共役の `Base`・`degFr` を手で計算する必要は無い)。

★★★**`𝒪^▷` の対応は独立の補題を要さない**。
`OTri P A = {φ | IsBaseIdentity P φ ∧ IsLinear P φ}` であり、
両方とも `Base` と `degFr` だけで決まるので、
`isFrobeniusNormalized_transport` の 3 仮定は**実質 2 つ**に落ちる。

★★仕上げは `pf_frobNormalizedType`(`𝒞^birat` 側、在庫)を流すだけ。
`IsOfFrobeniusNormalizedType (biratPre P G)` は
`IsOfBirationallyFrobeniusNormalizedType C P G` と**同じもの**である
(`BiratCat P G = C` かつ `(toBiratCat).obj A = A` なので)。

## ★☆見積りの推移(すべて「在庫を検索せずに測った」ことによる過大見積り)

**400 行 → 150 行 → 補題 2 本 → 実際は在庫の貼り合わせだけ**。
最後の 1 回は、`biratBase_biratPfHom` / `biratDeg_biratPfHom` の
「一般の元に対する版」を新規に作ろうとしたが、
★**`biratBase_biratPfHom_eq` / `biratDeg_biratPfHom_eq`(`Prop55PfSlice.lean`)が
まさにそれ**であった。教訓は毎回同じ: **在庫を検索してから見積もる。**
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

section PfBiratFrobNormalized

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} {P : PreFrobenioid C Φ} {F : FrobenioidCore P}
  [IsConnected D] {G : Frobenioid P}

omit [IsConnected D] in
/-- ★`𝒞` の対象がすべて isotropic なら `𝒞^birat` は Frobenius-isotropic 型。

★在庫の `birat_isOfFrobeniusIsotropicType`(`Prop55BiratOmega.lean`)は
`variable {D : Type u} [Category.{u2} D] {C : Type u2} [Category.{u2} C]` の節にあり
**宇宙が潰れている**ので、ここに一般宇宙版を置く(証明は同一)。 -/
theorem birat_isOfFrobeniusIsotropicType' (G : Frobenioid P) (hiso : ∀ W : C, IsIsotropic P W) :
    IsOfFrobeniusIsotropicType (biratPre P G) := fun A =>
  ⟨A, 𝟙 A, isFrobeniusType_of_isIso (biratPre P G) (𝟙 A),
    birat_isOfIsotropicType (G := G) hiso A⟩

/-! ## ★1. 第 1・2 段を束ねた `End` の乗法同型 -/

/-- ★★★★★**第 1・2 段を束ねた `End` の乗法同型** ——

    End_{(𝒞^birat)^pf}(A^{(n)})  ≃*  End_{(𝒞^pf)^birat}(Σ_{n·n} ⟨A^{(n)}, 1⟩)

★第 1 段は `biratPfHom`(`biratPfHom_comp` が乗法性)、
第 2 段は `scaleRootBirat`(関手なので `Functor.mapEnd` が `MonoidHom`)。
★★3 段目(同型による共役)は**束ねない** —— 在庫の
`isFrobeniusNormalized_of_iso` に任せた方が短い。 -/
noncomputable def pfBiratEndMulEquiv12 (hfi : IsOfFrobeniusIsotropicType P)
    (hiso : ∀ X : C, IsIsotropic P X)
    (Gpf : Frobenioid (pfRootPre P F)) (F' : FrobenioidCore (biratPre P G))
    (A : C) (n : ℕ+) :
    End (show PfCat (biratPre P G) F' from (biratUp P G (rtObj P F A n)))
      ≃* End ((scaleRootBirat (F := F) (n * n) Gpf).obj
        (show BiratCat (pfRootPre P F) Gpf from (⟨rtObj P F A n, 1⟩ : PfRootObj P F))) :=
  (biratPfEndMulEquiv hfi hiso Gpf F' (rtObj P F A n)).trans
    (scaleRootBiratEndMulEquiv (F := F) (n * n) Gpf
      (show BiratCat (pfRootPre P F) Gpf from (⟨rtObj P F A n, 1⟩ : PfRootObj P F)))

/-- ★★★★**束ねた乗法同型は底を保つ**(第 1・2 段の保存を継いだだけ)。 -/
theorem biratBase_pfBiratEndMulEquiv12 (hfi : IsOfFrobeniusIsotropicType P)
    (hiso : ∀ X : C, IsIsotropic P X)
    (Gpf : Frobenioid (pfRootPre P F)) (F' : FrobenioidCore (biratPre P G))
    (A : C) (n : ℕ+)
    (φ : End (show PfCat (biratPre P G) F' from (biratUp P G (rtObj P F A n)))) :
    biratBase (P := pfRootPre P F) (G := Gpf)
        ((pfBiratEndMulEquiv12 hfi hiso Gpf F' A n φ : End _) : _ ⟶ _)
      = pfBase (P := biratPre P G) (F := F') ((φ : End _) : _ ⟶ _) :=
  (biratBase_scaleRootBirat (F := F) (n * n) Gpf
      (biratPfHom hfi Gpf F' (rtObj P F A n) (rtObj P F A n) φ)).trans
    (biratBase_biratPfHom_eq hfi Gpf F' (rtObj P F A n) (rtObj P F A n) φ)

/-- ★★★★**束ねた乗法同型は Frobenius 次数を保つ**。 -/
theorem biratDeg_pfBiratEndMulEquiv12 (hfi : IsOfFrobeniusIsotropicType P)
    (hiso : ∀ X : C, IsIsotropic P X)
    (Gpf : Frobenioid (pfRootPre P F)) (F' : FrobenioidCore (biratPre P G))
    (A : C) (n : ℕ+)
    (φ : End (show PfCat (biratPre P G) F' from (biratUp P G (rtObj P F A n)))) :
    biratDeg (P := pfRootPre P F) (G := Gpf)
        ((pfBiratEndMulEquiv12 hfi hiso Gpf F' A n φ : End _) : _ ⟶ _)
      = pfDeg (P := biratPre P G) (F := F') ((φ : End _) : _ ⟶ _) :=
  (biratDeg_scaleRootBirat (F := F) (n * n) Gpf
      (biratPfHom hfi Gpf F' (rtObj P F A n) (rtObj P F A n) φ)).trans
    (biratDeg_biratPfHom_eq hfi Gpf F' (rtObj P F A n) (rtObj P F A n) φ)

/-! ## ★2. 移送 —— まず `Σ_{n·n} ⟨A^{(n)}, 1⟩` で示す -/

set_option maxHeartbeats 1000000 in
/-- ★★★★★**根を 1 に落とした対象での Frobenius-normalized 性**。

★出発点は `pf_frobNormalizedType`(`𝒞^birat` の `pf` は Frobenius-normalized 型)、
移送は `isFrobeniusNormalized_transport`。
★★`𝒪^▷` の対応(第 3 仮定)は `IsBaseIdentity ∧ IsLinear` なので、
`Base` と `degFr` の対応から `and_congr` で出る —— 独立の補題は要らない。 -/
theorem pfRoot_birat_frobNormalized_scaled (hfi : IsOfFrobeniusIsotropicType P)
    (hiso : ∀ X : C, IsIsotropic P X)
    (Gpf : Frobenioid (pfRootPre P F)) (F' : FrobenioidCore (biratPre P G))
    (hbfn : IsOfBirationallyFrobeniusNormalizedType C P G) (A : C) (n : ℕ+) :
    IsFrobeniusNormalized (biratPre (pfRootPre P F) Gpf)
      ((scaleRootBirat (F := F) (n * n) Gpf).obj
        (show BiratCat (pfRootPre P F) Gpf from (⟨rtObj P F A n, 1⟩ : PfRootObj P F))) := by
  have hb := biratBase_pfBiratEndMulEquiv12 hfi hiso Gpf F' A n
  have hd := biratDeg_pfBiratEndMulEquiv12 hfi hiso Gpf F' A n
  have h1 := hb 1
  have h2 := hd 1
  rw [map_one] at h1 h2
  have hbi : ∀ φ : End (show PfCat (biratPre P G) F' from (biratUp P G (rtObj P F A n))),
      IsBaseIdentity (pfPre (biratPre P G) F') ((φ : End _) : _ ⟶ _)
        ↔ IsBaseIdentity (biratPre (pfRootPre P F) Gpf)
          ((pfBiratEndMulEquiv12 hfi hiso Gpf F' A n φ : End _) : _ ⟶ _) :=
    fun φ => ⟨fun h => ((hb φ).trans h).trans h1.symm,
      fun h => (hb φ).symm.trans (h.trans h1)⟩
  refine isFrobeniusNormalized_transport (pfPre (biratPre P G) F')
    (biratPre (pfRootPre P F) Gpf) (pfBiratEndMulEquiv12 hfi hiso Gpf F' A n)
    hbi (fun φ => hd φ) (fun α => and_congr (hbi α) ?_)
    (pf_frobNormalizedType (F := F') (birat_isOfFrobeniusIsotropicType' G hiso)
      (fun X => hbfn X) (biratUp P G (rtObj P F A n)))
  exact ⟨fun h => (hd α).trans h, fun h => (hd α).symm.trans h⟩

/-! ## ★3. 一般の根へ —— 同型で移す -/

set_option maxHeartbeats 1000000 in
/-- ★★★★★★**[FrdI] Proposition 5.5, (iii)** ——
`𝒞` が birationally Frobenius-normalized 型なら **`𝒞^pf` もそう**。

原文 (FrdI p.105):
> if C is of standard (respectively, rationally standard) type, then so are Cun-tr, Crlf.

★★これが `Theorem 6.4, (i)` の 5 圏のうち `𝒞^pf` を止めていた**最後の条**である。

★仮定は原文の常備仮定(`hfi` / `hiso`)と `𝒞` の birat-Frobenius-normalized 性だけで、
**`hut`(unit-trivial)は要らない**——そこが在庫の `pfRoot_isOfModelType_of_arb` との違いである。 -/
theorem pfRoot_biratFrobNormalizedType (hfi : IsOfFrobeniusIsotropicType P)
    (hiso : ∀ X : C, IsIsotropic P X)
    (Gpf : Frobenioid (pfRootPre P F)) (F' : FrobenioidCore (biratPre P G))
    (hbfn : IsOfBirationallyFrobeniusNormalizedType C P G) :
    IsOfBirationallyFrobeniusNormalizedType (PfRootObj P F) (pfRootPre P F) Gpf := by
  rintro ⟨A, n⟩
  obtain ⟨ε, hε⟩ := pfRoot_exists_iso_root (F := F) A n n ((n * n) * 1) (by rw [mul_one])
  haveI := hε
  refine isFrobeniusNormalized_of_iso
    (((toBiratCat (pfRootPre P F) Gpf).mapIso (asIso ε)).symm) ?_
  exact pfRoot_birat_frobNormalized_scaled hfi hiso Gpf F' hbfn A n

end PfBiratFrobNormalized

/-! ### ★出典の紐付け -/

def pfBiratEndMulEquiv12.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 105,
    item := "Proposition 5.5, (iii) — End の乗法同型(第 1・2 段)",
    sectionId := "frdi-prop-5-5" }

def pfRoot_biratFrobNormalizedType.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 105,
    item := "Proposition 5.5, (iii) — 𝒞^pf は birationally Frobenius-normalized 型",
    sectionId := "frdi-prop-5-5" }

end ABC3.Found.FrdI
