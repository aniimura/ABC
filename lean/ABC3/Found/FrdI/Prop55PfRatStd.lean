/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Prop55PfBiratFn
import ABC3.Found.FrdI.Def24PfTransport

/-!
# [FrdI] Proposition 5.5, (iii) —— `𝒞^pf` は rationally standard 型

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.105。

原文 (FrdI p.105):
> if C is of standard (respectively, rationally standard) type, then so are Cun-tr, Crlf.

## ★★4 条をすべて集める

`Definition 4.5, (iii)` の 4 条が `𝒞 → 𝒞^pf` へ移ることを、それぞれの出どころで:

| 条 | 出どころ |
|---|---|
| (a) birationally Frobenius-normalized | ★`pfRoot_biratFrobNormalizedType`(`Prop55PfBiratFn.lean`) |
| (a) rational | ★`pfRoot_isOfStrictlyRationalType`(本ファイル) |
| (a) standard | `pfRoot_standardType`(在庫) |
| (b) `((𝒞^pf)^un-tr)^birat` の Frobenius-compact 対象 | 仮引数で受ける(模型ごとに与える) |

★(b) を仮引数で受けるのは在庫の `unTr_isOfRationallyStandardType` と同じ流儀である。
算術では `unTr_pf_biratCompact_of_baseId`(`Thm64RatStd.lean`)が与える。

## ★★★rational の条が短く済む理由

`IsStrictlyRational P G ι A` は

    ∀ p, ∃ a b, toGp a - toGp b ∈ Φ^birat(A) ∧ p ∈ Supp_ι(a) ∧ p ∉ Supp_ι(b)

であり、★**`ι` は自由な仮引数**である(`IsPerfFactorialWith` の 11 条は要らない)。
したがって `𝒞^pf` へ移すのに要るのは次の 3 つだけ:

1. **素点の対応** —— `primeToPf_bijective`(在庫、`MonoidPrime.lean`)。
2. **台の対応** —— `mem_suppElt_iotaPf`(`Def24PfTransport.lean`)。
3. **`Φ^birat` の対応** —— `phiBiratAt_pf_mem`(在庫、`Prop55PfSlice.lean`)。
   ★これは `Proposition 5.5, (iv)` の単系の同定の片側そのものである。

★根 `n` の違いは効かない —— `Φ^birat` は**底だけ**で決まり(`phiBiratOn_base`)、
`pfRootToElem` は `root` を見ないので `⟨A,n⟩` と `⟨A,1⟩` の底は同じである。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

section PfRatStd

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} {P : PreFrobenioid C Φ} {F : FrobenioidCore P}
  [IsConnected D] {G : Frobenioid P}

omit [IsConnected D] in
/-- ★★★**`𝒞^pf` では `Φ^birat` は根に依らない** —— 底が同じだから。

★`phiBiratOn_base` は `Φ^birat` が `𝒟` の対象(＝底)だけで決まることを言い、
`pfRootToElem` は `root` を見ない。 -/
theorem mem_phiBiratAt_pfRoot_root_indep (Gpf : Frobenioid (pfRootPre P F))
    (hisoPf : ∀ Y : PfRootObj P F, IsIsotropic (pfRootPre P F) Y) (X : PfRootObj P F)
    (y : Gp ((Φ.pfOn (phiSharp P)).val (((pfRootPre P F).toElem.obj X).base))) :
    y ∈ phiBiratAt (pfRootPre P F) Gpf (show BiratCat (pfRootPre P F) Gpf from X)
      ↔ y ∈ phiBiratAt (pfRootPre P F) Gpf
          (show BiratCat (pfRootPre P F) Gpf from (⟨X.obj, 1⟩ : PfRootObj P F)) := by
  rw [← phiBiratOn_base Gpf hisoPf X, ← phiBiratOn_base Gpf hisoPf (⟨X.obj, 1⟩ : PfRootObj P F)]
  exact Iff.rfl

/-- ★★★★★**`𝒞` が strictly rational 型なら `𝒞^pf` もそう**
(`Definition 4.5, (ii)` の `𝒞^pf` への伝播)。

★台の対応は `Def24PfTransport.lean` の `iotaPf` / `mem_suppElt_iotaPf`、
`Φ^birat` の対応は在庫の `phiBiratAt_pf_mem`。
素点は `primeToPf` で 1 対 1 に対応する(`primeToPf_bijective`)。 -/
theorem pfRoot_isOfStrictlyRationalType (hfi : IsOfFrobeniusIsotropicType P)
    (Gpf : Frobenioid (pfRootPre P F)) (F' : FrobenioidCore (biratPre P G))
    (hisoPf : ∀ Y : PfRootObj P F, IsIsotropic (pfRootPre P F) Y)
    (ι : ∀ Y : D, Prime (Φ.val Y) → Pf (Φ.val Y) → NNReal)
    (h : IsOfStrictlyRationalType C P G ι) :
    IsOfStrictlyRationalType (PfRootObj P F) (pfRootPre P F) Gpf
      (fun Y => iotaPf (P.divisorial Y) (ι Y)) := by
  intro X p
  obtain ⟨p', hp'⟩ := (primeToPf_bijective
    (isTorsionFreeNaive_of_divisorial (P.divisorial (P.toElem.obj X.obj).base))).2 p
  subst hp'
  obtain ⟨a, b, hab, ha, hb⟩ := h X.obj p'
  refine ⟨Pf.of a, Pf.of b, ?_, ?_, fun hc => hb ?_⟩
  · refine (mem_phiBiratAt_pfRoot_root_indep Gpf hisoPf X _).mpr ?_
    have key := phiBiratAt_pf_mem (F := F) hfi Gpf F' X.obj hab
    rw [map_sub, gpMap_toGp, gpMap_toGp] at key
    exact key
  · exact (mem_suppElt_iotaPf _ (ι _) a p').mpr ha
  · exact (mem_suppElt_iotaPf _ (ι _) b p').mp hc

/-- ★★★★★★**[FrdI] Proposition 5.5, (iii)** —— `𝒞^pf` は **rationally standard 型**。

原文 (FrdI p.105):
> if C is of standard (respectively, rationally standard) type, then so are Cun-tr, Crlf.

★★4 条のうち 3 条は上の表のとおり在庫と本ファイル・`Prop55PfBiratFn.lean` で埋まる。
残る (b)(`((𝒞^pf)^un-tr)^birat` の Frobenius-compact 対象)は
在庫の `unTr_isOfRationallyStandardType` と同じく**仮引数で受ける** ——
模型ごとに与えるものであり(算術では `unTr_pf_biratCompact_of_baseId`)、
`𝒞` 側の rationally standard 性からは出ない。 -/
theorem pfRoot_isOfRationallyStandardType (hfi : IsOfFrobeniusIsotropicType P)
    (hiso : ∀ X : C, IsIsotropic P X)
    (hisoPf : ∀ Y : PfRootObj P F, IsIsotropic (pfRootPre P F) Y)
    (Gpf : Frobenioid (pfRootPre P F)) (F' : FrobenioidCore (biratPre P G))
    (ι : ∀ Y : D, Prime (Φ.val Y) → Pf (Φ.val Y) → NNReal)
    (hbfn : IsOfBirationallyFrobeniusNormalizedType C P G)
    (hsr : IsOfStrictlyRationalType C P G ι)
    (hngl : ¬ IsOfGroupLikeType P)
    (hstd : IsOfStandardType D C P F)
    (hcpt : ∃ X : BiratCat (unTrPre (pfRootPre P F) Gpf.core)
        (unTr_frobenioid (pfRootPre P F) Gpf.core Gpf),
      IsFrobeniusCompact (unTrBiratPre (pfRootPre P F) Gpf.core Gpf) X) :
    IsOfRationallyStandardType (pfRootPre P F) Gpf
      (fun Y => iotaPf (P.divisorial Y) (ι Y)) where
  biratFrobNormalized := pfRoot_biratFrobNormalizedType hfi hiso Gpf F' hbfn
  rational := fun X => isRational_of_isStrictlyRational _ X
    (pfRoot_isOfStrictlyRationalType hfi Gpf F' hisoPf ι hsr X)
  standard := pfRoot_standardType hfi Gpf.core hngl hstd
  unTrBiratCompact := hcpt

end PfRatStd

/-! ### ★出典の紐付け -/

def pfRoot_isOfStrictlyRationalType.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 105,
    item := "Proposition 5.5, (iii) — 𝒞^pf は strictly rational 型",
    sectionId := "frdi-prop-5-5" }

def pfRoot_isOfRationallyStandardType.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 105,
    item := "Proposition 5.5, (iii) — 𝒞^pf は rationally standard 型",
    sectionId := "frdi-prop-5-5" }

def pfRoot_isOfRationallyStandardType.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "pfRoot_biratFrobNormalizedType((a) birat-Frobenius-normalized)"
      (.inProject "ABC3" "ABC3.Found.FrdI.pfRoot_biratFrobNormalizedType") 105,
    .citation "[ABC3]" "pfRoot_standardType((a) standard)"
      (.inProject "ABC3" "ABC3.Found.FrdI.pfRoot_standardType") 105,
    .citation "[ABC3]" "phiBiratAt_pf_mem(Φ^birat の像は (Φ^pf)^birat に入る)"
      (.inProject "ABC3" "ABC3.Found.FrdI.phiBiratAt_pf_mem") 105,
    .citation "[ABC3]" "mem_suppElt_iotaPf / primeToPf_bijective(台と素点の対応)"
      (.inProject "ABC3" "ABC3.Found.FrdI.mem_suppElt_iotaPf") 105,
    .implicitStep
      ("★(b)(`((𝒞^pf)^un-tr)^birat` の Frobenius-compact 対象)は仮引数で受ける。" ++
       "在庫の unTr_isOfRationallyStandardType と同じ流儀で、模型ごとに与えるもの" ++
       "(算術では unTr_pf_biratCompact_of_baseId)") 105 ]

end ABC3.Found.FrdI
