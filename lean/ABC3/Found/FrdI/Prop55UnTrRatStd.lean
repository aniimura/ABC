/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Prop55PfRatStd
import ABC3.Found.FrdI.Thm64RatStd
import ABC3.Found.FrdI.Prop53UntrBirat
import ABC3.Found.FrdI.Thm52ModelType

/-!
# [FrdI] Proposition 5.5, (iii) —— `𝒞^un-tr` と `(𝒞^pf)^un-tr` は rationally standard 型

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.105。

原文 (FrdI p.105):
> if C is of standard (respectively, rationally standard) type, then so are Cun-tr, Crlf.

## ★★在庫の `unTr_isOfRationallyStandardType` との違い

在庫のものは `ModelData` を要求する(`h : M.Hyp`)ので **`𝒞^pf` には当たらない** ——
`𝒞^pf` が model 型であることは `hut`(unit-trivial)を要し、
`Example 6.1` / `Example 6.3` の `𝒞` はそれを満たさないからである。

★本ファイルは同じ結論を**一般形**で出す。4 条の出どころは:

| 条 | 出どころ | 一般性 |
|---|---|---|
| (a) birat-Frobenius-normalized | `unTr_isOfModelType`(在庫) | ★一般 |
| (a) rational | ★`unTr_isOfStrictlyRationalType`(本ファイル) | ★一般 |
| (a) standard | `unTr_standardType`(在庫) | ★一般 |
| (b) Frobenius-compact 対象 | 仮引数 | 模型ごと |

★rational が一般に出るのは **`Φ^birat` が `𝒞` と `𝒞^un-tr` で点ごとに一致する**
(`phiBiratAt_unTr_eq`、在庫)からで、`Φ` も `ι` もそのまま使える。

## ★★★`(𝒞^pf)^un-tr` —— `Theorem 6.4, (i)` の 5 圏のうち 4 つ目

上の一般形を **`𝒞^pf` に当てる**だけである。(b) の対象は
`(((𝒞^pf)^un-tr)^un-tr)^birat` の中に要るので、`Φ^birat` を 3 段たどる:

    Φ^birat(((𝒞^pf)^un-tr)^un-tr)  =  Φ^birat((𝒞^pf)^un-tr)   (phiBiratAt_unTr_eq)
                                   =  Φ^birat(𝒞^pf)           (phiBiratAt_unTr_eq)
                                   ∋  Pf.of(Φ^birat(𝒞) の元)  (phiBiratAt_pf_mem)

★底の自己射が恒等である対象(算術では `⊥` = `F` 自身)を取れば
`Φ^pf` の側でも恒等になる(`Pf.map (id) = id`)ので、`hbase` はそのまま持ち上がる。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

section UnTrRatStd

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} {P : PreFrobenioid C Φ} [IsConnected D]

omit [IsConnected D] in
/-- ★★★★★**`𝒞` が strictly rational 型なら `𝒞^un-tr` もそう**(★一般形)。

★`Φ` は `𝒞` と `𝒞^un-tr` で**同じ**なので `ι` はそのまま、
`Φ^birat` も点ごとに一致する(`phiBiratAt_unTr_eq`)。
★★在庫の `unTr_isOfStrictlyRationalType_of_hsp` は `ModelData` 経由だったが、
中身は `phiBiratAt_unTr_eq` 1 本なので model は要らない。 -/
theorem unTr_isOfStrictlyRationalType (Fc : FrobenioidCore P) (G : Frobenioid P)
    (hiso : ∀ X : C, IsIsotropic P X)
    (ι : ∀ Y : D, Prime (Φ.val Y) → Pf (Φ.val Y) → NNReal)
    (h : IsOfStrictlyRationalType C P G ι) :
    IsOfStrictlyRationalType (UnTr P) (unTrPre P Fc) (unTr_frobenioid P Fc G) ι := by
  intro X p
  obtain ⟨a, b, hmem, ha, hb⟩ := h X.obj p
  refine ⟨a, b, ?_, ha, hb⟩
  rw [phiBiratAt_unTr_eq Fc G hiso (show Istr P from X)]
  exact hmem

/-- ★★★★★★**[FrdI] Proposition 5.5, (iii)** —— `𝒞^un-tr` は rationally standard 型
(★**`ModelData` を要求しない一般形**)。

★(b)(`((𝒞^un-tr)^un-tr)^birat` の Frobenius-compact 対象)だけは仮引数で受ける ——
模型ごとに与えるものであり、`𝒞` 側の rationally standard 性からは出ない。 -/
theorem unTr_isOfRationallyStandardType_gen (Fc : FrobenioidCore P) (G : Frobenioid P)
    (hiso : ∀ X : C, IsIsotropic P X)
    (ι : ∀ Y : D, Prime (Φ.val Y) → Pf (Φ.val Y) → NNReal)
    (hsr : IsOfStrictlyRationalType C P G ι)
    (hint : ∀ A : D, IsIntegralMonoid (Φ.val A))
    (hfsmD : IsOfFSMType D)
    (hngl : ¬ IsOfGroupLikeType P)
    (hstd : IsOfStandardType D C P Fc)
    (hcpt : ∃ X : BiratCat (unTrPre (unTrPre P Fc) (unTr_frobenioid P Fc G).core)
        (unTr_frobenioid (unTrPre P Fc) (unTr_frobenioid P Fc G).core
          (unTr_frobenioid P Fc G)),
      IsFrobeniusCompact (unTrBiratPre (unTrPre P Fc) (unTr_frobenioid P Fc G).core
        (unTr_frobenioid P Fc G)) X) :
    IsOfRationallyStandardType (unTrPre P Fc) (unTr_frobenioid P Fc G) ι where
  biratFrobNormalized := (unTr_isOfModelType Fc G).2
  rational := fun X => isRational_of_isStrictlyRational ι X
    (unTr_isOfStrictlyRationalType Fc G hiso ι hsr X)
  standard := unTr_standardType Fc G hint hfsmD _ hngl hstd
  unTrBiratCompact := hcpt

end UnTrRatStd

section PfUnTr

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} {P : PreFrobenioid C Φ} {F : FrobenioidCore P}
  [IsConnected D] {G : Frobenioid P}

set_option maxHeartbeats 1000000 in
/-- ★★★★★**`(((𝒞^pf)^un-tr)^un-tr)^birat` の Frobenius-compact 対象**。

★`Φ^birat` を 3 段たどるだけである(ファイル冒頭の式)。
底の自己射が恒等な対象を `𝒞` で取れば、`Φ^pf` の側でも恒等になる
(`Pf.map (AddMonoidHom.id _) = AddMonoidHom.id _`)。 -/
theorem pf_unTr2_biratCompact_of_baseId (hfi : IsOfFrobeniusIsotropicType P)
    (hisoPf : ∀ Y : PfRootObj P F, IsIsotropic (pfRootPre P F) Y)
    (Gpf : Frobenioid (pfRootPre P F)) (F' : FrobenioidCore (biratPre P G)) (A : C)
    (hb : ∀ e : (P.toElem.obj A).base ⟶ (P.toElem.obj A).base,
      Φ.map e = AddMonoidHom.id (Φ.val (P.toElem.obj A).base))
    (z : Gp (Φ.val (P.toElem.obj A).base))
    (hzmem : z ∈ phiBiratAt P G (biratUp P G A))
    (hz : ∀ n : ℕ, 0 < n →
      n • gpMap (Φ.val (P.toElem.obj A).base) (Pf.of) z ≠ 0) :
    ∃ X : BiratCat
        (unTrPre (unTrPre (pfRootPre P F) Gpf.core)
          (unTr_frobenioid (pfRootPre P F) Gpf.core Gpf).core)
        (unTr_frobenioid (unTrPre (pfRootPre P F) Gpf.core)
          (unTr_frobenioid (pfRootPre P F) Gpf.core Gpf).core
          (unTr_frobenioid (pfRootPre P F) Gpf.core Gpf)),
      IsFrobeniusCompact (unTrBiratPre (unTrPre (pfRootPre P F) Gpf.core)
        (unTr_frobenioid (pfRootPre P F) Gpf.core Gpf).core
        (unTr_frobenioid (pfRootPre P F) Gpf.core Gpf)) X := by
  set Q := unTrPre (pfRootPre P F) Gpf.core with hQ
  set GQ := unTr_frobenioid (pfRootPre P F) Gpf.core Gpf with hGQ
  set X1 : Istr (pfRootPre P F) := ⟨(⟨A, 1⟩ : PfRootObj P F), hisoPf _⟩ with hX1
  set X2 : Istr Q := ⟨show UnTr (pfRootPre P F) from X1,
    unTr_isotropic (pfRootPre P F) Gpf.core _⟩ with hX2
  have hpf : ∀ e : (P.toElem.obj A).base ⟶ (P.toElem.obj A).base,
      (Φ.pfOn (phiSharp P)).map e
        = AddMonoidHom.id ((Φ.pfOn (phiSharp P)).val (P.toElem.obj A).base) := by
    intro e
    show Pf.map (Φ.map e) = _
    rw [hb e]
    exact Pf.map_id
  have e1 : phiBiratAt (unTrPre Q GQ.core) (unTr_frobenioid Q GQ.core GQ)
      (show UnTr Q from X2) = phiBiratAt Q GQ X2.obj :=
    phiBiratAt_unTr_eq GQ.core GQ (fun Y => unTr_isotropic (pfRootPre P F) Gpf.core Y) X2
  have e2 : phiBiratAt Q GQ (show UnTr (pfRootPre P F) from X1)
      = phiBiratAt (pfRootPre P F) Gpf X1.obj :=
    phiBiratAt_unTr_eq Gpf.core Gpf hisoPf X1
  have e3 := phiBiratAt_pf_mem (F := F) hfi Gpf F' A hzmem
  refine ⟨show UnTr Q from X2, birat_isFrobeniusCompact_of_baseId (unTrPre Q GQ.core)
    (unTr_frobenioid Q GQ.core GQ)
    (fun Z => (unTr_isOfModelType GQ.core GQ).2 Z)
    (show UnTr Q from X2) (fun θ => hpf _)
    (gpMap (Φ.val (P.toElem.obj A).base) (Pf.of) z) ?_ hz⟩
  rw [e1]
  show _ ∈ phiBiratAt Q GQ (show UnTr (pfRootPre P F) from X1)
  rw [e2]
  exact e3

set_option maxHeartbeats 1000000 in
/-- ★★★★★★**[FrdI] Proposition 5.5, (iii)** —— **`(𝒞^pf)^un-tr` は rationally standard 型**。

★★これで `Theorem 6.4, (i)` が並べる 5 圏のうち **4 つ**が揃う
(`𝒞`・`𝒞^un-tr`・`𝒞^pf`・`(𝒞^pf)^un-tr`)。残るは `𝒞^rlf` である。

★原文は `(𝒞^pf)^un-tr` を `Proposition 5.5, (iii)` に委ねているが、
在庫の `unTr_isOfRationallyStandardType` は `ModelData` を要求するので当たらない
(`𝒞^pf` の model 型は `hut` を要する)。★本定理は一般形で回り込む。 -/
theorem pfRoot_unTr_isOfRationallyStandardType (hfi : IsOfFrobeniusIsotropicType P)
    (hisoPf : ∀ Y : PfRootObj P F, IsIsotropic (pfRootPre P F) Y)
    (Gpf : Frobenioid (pfRootPre P F)) (F' : FrobenioidCore (biratPre P G))
    (ι : ∀ Y : D, Prime (Φ.val Y) → Pf (Φ.val Y) → NNReal)
    (hsr : IsOfStrictlyRationalType C P G ι)
    (hfsmD : IsOfFSMType D)
    (hngl : ¬ IsOfGroupLikeType P)
    (hstd : IsOfStandardType D C P F)
    (A : C)
    (hb : ∀ e : (P.toElem.obj A).base ⟶ (P.toElem.obj A).base,
      Φ.map e = AddMonoidHom.id (Φ.val (P.toElem.obj A).base))
    (z : Gp (Φ.val (P.toElem.obj A).base))
    (hzmem : z ∈ phiBiratAt P G (biratUp P G A))
    (hz : ∀ n : ℕ, 0 < n → n • gpMap (Φ.val (P.toElem.obj A).base) (Pf.of) z ≠ 0) :
    IsOfRationallyStandardType (unTrPre (pfRootPre P F) Gpf.core)
      (unTr_frobenioid (pfRootPre P F) Gpf.core Gpf)
      (fun Y => iotaPf (P.divisorial Y) (ι Y)) :=
  unTr_isOfRationallyStandardType_gen Gpf.core Gpf hisoPf _
    (pfRoot_isOfStrictlyRationalType hfi Gpf F' hisoPf ι hsr)
    (fun Y => (Pf.isDivisorial' (P.divisorial Y)).1.1)
    hfsmD
    (pfRoot_not_isOfGroupLikeType hngl)
    (pfRoot_standardType hfi Gpf.core hngl hstd)
    (pf_unTr2_biratCompact_of_baseId hfi hisoPf Gpf F' A hb z hzmem hz)

end PfUnTr

/-! ### ★出典の紐付け -/

def unTr_isOfStrictlyRationalType.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 105,
    item := "Proposition 5.5, (iii) — 𝒞^un-tr は strictly rational 型(一般形)",
    sectionId := "frdi-prop-5-5" }

def unTr_isOfRationallyStandardType_gen.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 105,
    item := "Proposition 5.5, (iii) — 𝒞^un-tr は rationally standard 型(一般形)",
    sectionId := "frdi-prop-5-5" }

def pfRoot_unTr_isOfRationallyStandardType.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 105,
    item := "Proposition 5.5, (iii) — (𝒞^pf)^un-tr は rationally standard 型",
    sectionId := "frdi-prop-5-5" }

def pfRoot_unTr_isOfRationallyStandardType.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "phiBiratAt_unTr_eq(Φ^birat は 𝒞 と 𝒞^un-tr で一致)"
      (.inProject "ABC3" "ABC3.Found.FrdI.phiBiratAt_unTr_eq") 105,
    .citation "[ABC3]" "phiBiratAt_pf_mem(Φ^birat の像は (Φ^pf)^birat に入る)"
      (.inProject "ABC3" "ABC3.Found.FrdI.phiBiratAt_pf_mem") 105,
    .citation "[ABC3]" "unTr_isOfModelType / unTr_standardType(在庫。どちらも一般)"
      (.inProject "ABC3" "ABC3.Found.FrdI.unTr_isOfModelType") 105,
    .citation "[ABC3]" "pfRoot_isOfStrictlyRationalType / pfRoot_standardType(𝒞^pf の側)"
      (.inProject "ABC3" "ABC3.Found.FrdI.pfRoot_isOfStrictlyRationalType") 105,
    .implicitStep
      ("★在庫の unTr_isOfRationallyStandardType は ModelData を要求するので " ++
       "𝒞^pf には当たらない(𝒞^pf の model 型は hut を要する)。一般形で回り込んだ") 105 ]

end ABC3.Found.FrdI
