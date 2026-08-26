/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Prop53Birat
import ABC3.Found.FrdI.Prop48Cond3
import ABC3.Found.FrdI.Prop48Cpt
import ABC3.Found.FrdI.Thm52NatIso
import ABC3.Found.FrdI.Thm52Frob

/-!
# [FrdI] Theorem 6.2, (iii) —— rationally standard 型の組み立て

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.111。

原文 (FrdI p.111):
> nioid C is of isotropic, standard, and birationally Frobenius-normalized

原文 (FrdI p.111):
> of an element of B(L), then C is of rationally standard type.

## ★★何が残っていたか

`Definition 4.5, (iii)` の `IsOfRationallyStandardType` は 4 条:

| 条 | 状況 |
|---|---|
| (a) birationally Frobenius-normalized 型 | ★model では**自動**(`model_isOfBiratFrobNormalizedType`) |
| (a) rational 型 | `isOfRationalType_of_divB`(`Prop53Birat.lean`、済) |
| (a) standard 型 | `model_isOfStandardType`(`Thm52Frob.lean`、6 条中 3 条は自動) |
| (b) `(𝒞^un-tr)^birat` の Frobenius-compact 対象 | ★**本ファイル** |

★(b) の中身は `birat_isFrobeniusCompact_of_ne_zero`(`Prop48Cond3.lean`、在庫)であり、
入力は「`Div : 𝒪^▷ ↠ Φ`」「`(𝒞^un-tr)^birat` が Frobenius-normalized」
「`Φ` が non-dilating」「`Φ` に 0 でない元が 1 つある」の 4 つである。

★★**`Div : 𝒪^▷ ↠ Φ` は `𝒞^istr` の側で持っていれば足りる**
(`unTr_divSurj_of_istr`)—— `𝒞^un-tr` は射を `𝔽_Φ` の像で割った商なので、
`Base` も `degFr` も `Div` も**代表元のものがそのまま降りる**(3 つとも `rfl` で済んだ)。

★★「`(𝒞^un-tr)^birat` が Frobenius-normalized」は
`unTr_isOfModelType` の第 2 成分がそのまま与える(在庫)。

## ★本ファイルで閉じること

| 定理 | 中身 |
|---|---|
| `unTr_divSurj_of_istr` | ★`Div : 𝒪^▷ ↠ Φ` が `𝒞^istr → 𝒞^un-tr` で降りる |
| `unTrBiratCompact_of_ne_zero` | ★★`(𝒞^un-tr)^birat` の Frobenius-compact 対象 |
| `ModelData.model_isOfRationallyStandardType` | ★★★★★★**model Frobenioid は rationally standard 型** |
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} (P : PreFrobenioid C Φ) (Fc : FrobenioidCore P)

/-! ## ★1. `Div : 𝒪^▷ ↠ Φ` は `𝒞^un-tr` へ降りる -/

/-- ★★**`Div : 𝒪^▷ ↠ Φ` は `𝒞^istr` から `𝒞^un-tr` へ移る**。

★`𝒞^un-tr` は射を `𝔽_Φ` の像で割った商なので、
`Base`・`degFr`・`Div` はどれも**代表元のものがそのまま降りる**。 -/
theorem unTr_divSurj_of_istr
    (hdivS : ∀ (Y : Istr P) (a : Φ.val ((istrPre P Fc).toElem.obj Y).base),
      ∃ u : OTri (istrPre P Fc) Y, (istrPre P Fc).Div ((u : End Y) : Y ⟶ Y) = a)
    (Y : UnTr P) (a : Φ.val ((unTrPre P Fc).toElem.obj Y).base) :
    ∃ u : OTri (unTrPre P Fc) Y, (unTrPre P Fc).Div ((u : End Y) : Y ⟶ Y) = a := by
  obtain ⟨u, hu⟩ := hdivS (show Istr P from Y) a
  refine ⟨⟨(istrToUnTr P).map ((u : End (show Istr P from Y)) : _ ⟶ _), ?_, ?_⟩, hu⟩
  · show (unTrPre P Fc).Base ((istrToUnTr P).map ((u : End (show Istr P from Y)) : _ ⟶ _))
      = (unTrPre P Fc).Base (𝟙 Y)
    exact u.2.1
  · show (unTrPre P Fc).degFr ((istrToUnTr P).map ((u : End (show Istr P from Y)) : _ ⟶ _)) = 1
    exact u.2.2

/-! ## ★2. `(𝒞^un-tr)^birat` の Frobenius-compact 対象 -/

/-- ★★★**`Definition 4.5, (iii), (b)`** —— `(𝒞^un-tr)^birat` の Frobenius-compact 対象。

★中身は `birat_isFrobeniusCompact_of_ne_zero`(在庫)を `𝒞^un-tr` に当てるだけ。 -/
theorem unTrBiratCompact_of_ne_zero (G : Frobenioid P)
    (hdivS : ∀ (Y : Istr P) (a : Φ.val ((istrPre P Fc).toElem.obj Y).base),
      ∃ u : OTri (istrPre P Fc) Y, (istrPre P Fc).Div ((u : End Y) : Y ⟶ Y) = a)
    (hfn : ∀ X : BiratCat (unTrPre P Fc) (unTr_frobenioid P Fc G),
      IsFrobeniusNormalized (unTrBiratPre P Fc G) X)
    (hndOn : MonoidOn.IsNonDilatingOn Φ)
    (A : BiratCat (unTrPre P Fc) (unTr_frobenioid P Fc G))
    (x₀ : Φ.val ((unTrPre P Fc).toElem.obj
      (biratDown (unTrPre P Fc) (unTr_frobenioid P Fc G) A)).base) (hx₀ : x₀ ≠ 0) :
    ∃ X : BiratCat (unTrPre P Fc) (unTr_frobenioid P Fc G),
      IsFrobeniusCompact (unTrBiratPre P Fc G) X :=
  ⟨A, birat_isFrobeniusCompact_of_ne_zero (unTr_divSurj_of_istr P Fc hdivS) hfn hndOn A x₀ hx₀⟩

/-! ## ★3. model Frobenioid は rationally standard 型 -/

namespace ModelData

variable {M : ModelData.{v, u, w} D}

/-- ★★★★★★**[FrdI] Theorem 6.2, (iii)** —— model Frobenioid が **rationally standard 型**。

★4 条のうち **(a) の第 1 条は model では自動**である
(`model_isOfBiratFrobNormalizedType`、`model_ratStandardType_iff` が吸収済み)。

原文 (FrdI p.111):
> of an element of B(L), then C is of rationally standard type. -/
theorem model_isOfRationallyStandardType (h : Hyp M)
    (ι : ∀ Y : D, Prime (M.phi.val Y) → Pf (M.phi.val Y) → NNReal)
    (R : RatFnData (modelPre h) (model_frobenioid h))
    (hsp : ∀ (A : Obj M) (p : Prime (M.phi.val ((modelPre h).toElem.obj A).base)),
      ∃ (a b : M.phi.val ((modelPre h).toElem.obj A).base)
        (y : R.bmon.val ((modelPre h).toElem.obj A).base),
        (toGp _ a - toGp _ b = R.divB _ y ∨ toGp _ a - toGp _ b = -(R.divB _ y)) ∧
        p ∈ SuppElt (ι _) a ∧ p ∉ SuppElt (ι _) b)
    (hfsmff : IsOfFSMFFType D) (hnd : M.phi.IsNonDilatingOn)
    (hgl : IsOfGroupLikeType (modelPre h) →
      ∃ A, IsFrobeniusCompact (istrPre (modelPre h) (model_frobenioid h).core) A)
    (hdivS : ∀ (Y : Istr (modelPre h))
        (a : M.phi.val ((istrPre (modelPre h) (model_frobenioid h).core).toElem.obj Y).base),
      ∃ u : OTri (istrPre (modelPre h) (model_frobenioid h).core) Y,
        (istrPre (modelPre h) (model_frobenioid h).core).Div ((u : End Y) : Y ⟶ Y) = a)
    (A : BiratCat (unTrPre (modelPre h) (model_frobenioid h).core)
      (unTr_frobenioid (modelPre h) (model_frobenioid h).core (model_frobenioid h)))
    (x₀ : M.phi.val ((unTrPre (modelPre h) (model_frobenioid h).core).toElem.obj
      (biratDown _ _ A)).base) (hx₀ : x₀ ≠ 0) :
    haveI := h.connectedD
    IsOfRationallyStandardType (modelPre h) (model_frobenioid h) ι :=
  haveI := h.connectedD
  (model_ratStandardType_iff h ι).mpr
    ⟨isOfRationalType_of_divB R ι hsp,
     model_isOfStandardType h (model_frobenioid h).core hfsmff hnd hgl,
     unTrBiratCompact_of_ne_zero (modelPre h) (model_frobenioid h).core (model_frobenioid h)
       hdivS
       (fun X => (unTr_isOfModelType (model_frobenioid h).core (model_frobenioid h)).2 X)
       hnd A x₀ hx₀⟩

end ModelData

/-! ### ★出典の紐付け -/

/-- ★★★★★locator —— `Theorem 6.2, (iii)` の rationally standard 型。 -/
def ModelData.model_isOfRationallyStandardType.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 111,
    item := "Theorem 6.2, (iii) — model Frobenioid は rationally standard 型",
    sectionId := "frdi-thm-6-2" }

/-- ★★★locator —— `Definition 4.5, (iii), (b)` の Frobenius-compact 対象。 -/
def unTrBiratCompact_of_ne_zero.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 111,
    item := "Theorem 6.2, (iii) — (𝒞^un-tr)^birat の Frobenius-compact 対象",
    sectionId := "frdi-thm-6-2" }

/-! ## ★4. 逸脱 (B) を外す —— 「主因子」ではなく「台に入る」で足りる

★★在庫の `birat_isFrobeniusCompact_of_ne_zero` は
**逸脱 (B)**(`hdivS : Div : 𝒪^▷(A) ↠ Φ(A)`)に乗っている。
★★★**この仮定は `Example 6.1` では成り立たない** —— 模型 Frobenioid では
`𝒪^▷(A)`(底が恒等・線型な自己射)は `cond` から `Div(φ) = Div_B(u_φ)` を満たすので、
`Div(𝒪^▷(A))` は **`Φ(A) ∩ Div_B(B(A))`(＝主因子)**であって `Φ(A)` 全体ではない。

★`hdivS` が使われるのは **1 箇所だけ**である ——
`phiBiratAt_eq_top_of_divSurj` を経由して「`toGp x ∈ phiBiratAt`」を得るところ。
★★したがって**必要なのは「その元が `Φ^birat` に入る」だけ**であり、
`birat_cd_eq` の側は `addMonoidHom_eq_id_of_primary_mprec` に差し替えれば
★**準素元についてだけ要求すれば足りる**。

★これが原文 `Theorem 6.2, (iii)` の追加仮定
「どの `D ∈ D_L` も `B(L)` の元の像の**台**に入る」の正しい使い方である。 -/

section NoDeviation

variable (G : Frobenioid P)

/-- ★★★★**`c = d`(逸脱 (B) を外した版)** —— `phiBiratAt` に入るのは**準素元だけでよい**。 -/
theorem birat_cd_eq_of_primary
    (A : BiratCat P G) (θ : A ≅ A) (c d : ℕ+)
    (hnd : IsNonDilating (Φ.map (biratBase θ.inv)))
    (hprim : ∀ x : Φ.val (P.toElem.obj (biratDown P G A)).base, IsPrimaryElt x →
      toGp _ x ∈ phiBiratAt P G A)
    (x₀ : Φ.val (P.toElem.obj (biratDown P G A)).base)
    (hx₀mem : toGp _ x₀ ∈ phiBiratAt P G A)
    (hx₀ : ∀ n : ℕ, 0 < n →
      n • toGp (Φ.val (P.toElem.obj (biratDown P G A)).base) x₀ ≠ 0)
    (hyp : ∀ u : End A, u ∈ OTimes (biratPre P G) A → ∃ k : ℕ+,
      ((endConj θ u) ^ (((d : ℕ+) : ℕ) * ((k : ℕ+) : ℕ)) : End A)
        = (u ^ (((c : ℕ+) : ℕ) * ((k : ℕ+) : ℕ)) : End A)) :
    c = d := by
  classical
  have hkey : ∀ x : Φ.val (P.toElem.obj (biratDown P G A)).base,
      toGp _ x ∈ phiBiratAt P G A →
      ∃ k : ℕ+, (((d : ℕ+) : ℕ) * ((k : ℕ+) : ℕ))
          • gpMap _ (Φ.map (biratBase θ.inv))
              (toGp (Φ.val (P.toElem.obj (biratDown P G A)).base) x)
        = (((c : ℕ+) : ℕ) * ((k : ℕ+) : ℕ))
          • toGp (Φ.val (P.toElem.obj (biratDown P G A)).base) x := by
    intro x hmem
    obtain ⟨u, hu, hux⟩ := hmem
    obtain ⟨k, hk⟩ := hyp u hu
    refine ⟨k, ?_⟩
    have hcj : (endConj θ u) ∈ OTimes (biratPre P G) A :=
      endConj_mem_otimes (biratPre P G) θ hu
    have h := congrArg (fun t : End A => biratDivGp ((t : A ⟶ A))) hk
    rw [biratDivGp_pow hcj, biratDivGp_pow hu, biratDivGp_endConj θ hu, hux] at h
    exact h
  have hsharp : IsSharp (Φ.val (P.toElem.obj (biratDown P G A)).base) :=
    (P.divisorial _).2
  have hint : IsIntegralMonoid (Φ.val (P.toElem.obj (biratDown P G A)).base) :=
    (P.divisorial _).1.1
  have hσ : Φ.map (biratBase θ.inv) = AddMonoidHom.id _ := by
    refine addMonoidHom_eq_id_of_primary_mprec hsharp _ hnd (fun x hpx => ?_)
    obtain ⟨k, hk⟩ := hkey x (hprim x hpx)
    have hmono : (((d : ℕ+) : ℕ) * ((k : ℕ+) : ℕ)) • (Φ.map (biratBase θ.inv) x)
        = (((c : ℕ+) : ℕ) * ((k : ℕ+) : ℕ)) • x := by
      refine hint ?_
      rw [toGp_nsmul, toGp_nsmul, ← gpMap_toGp _ (Φ.map (biratBase θ.inv)) x]
      exact hk
    refine ⟨((c : ℕ+) : ℕ) * ((k : ℕ+) : ℕ), Nat.mul_pos c.2 k.2, ?_⟩
    refine ⟨(((d : ℕ+) : ℕ) * ((k : ℕ+) : ℕ) - 1) • (Φ.map (biratBase θ.inv) x), ?_⟩
    have hpos : 0 < ((d : ℕ+) : ℕ) * ((k : ℕ+) : ℕ) := Nat.mul_pos d.2 k.2
    have h1 : Φ.map (biratBase θ.inv) x
          + (((d : ℕ+) : ℕ) * ((k : ℕ+) : ℕ) - 1) • (Φ.map (biratBase θ.inv) x)
        = (1 + (((d : ℕ+) : ℕ) * ((k : ℕ+) : ℕ) - 1)) • (Φ.map (biratBase θ.inv) x) := by
      rw [add_nsmul, one_nsmul]
    rw [h1, show 1 + (((d : ℕ+) : ℕ) * ((k : ℕ+) : ℕ) - 1)
        = ((d : ℕ+) : ℕ) * ((k : ℕ+) : ℕ) from by omega]
    exact hmono
  obtain ⟨k, hk⟩ := hkey x₀ hx₀mem
  rw [hσ] at hk
  have hk' : (((d : ℕ+) : ℕ) * ((k : ℕ+) : ℕ))
        • toGp (Φ.val (P.toElem.obj (biratDown P G A)).base) x₀
      = (((c : ℕ+) : ℕ) * ((k : ℕ+) : ℕ))
        • toGp (Φ.val (P.toElem.obj (biratDown P G A)).base) x₀ := by
    have hid : gpMap _ (AddMonoidHom.id (Φ.val (P.toElem.obj (biratDown P G A)).base))
        = AddMonoidHom.id _ := gpMap_id _
    rw [hid] at hk
    exact hk
  by_contra hne
  have hkpos : 0 < ((k : ℕ+) : ℕ) := k.2
  rcases Nat.lt_or_ge (((d : ℕ+) : ℕ) * ((k : ℕ+) : ℕ))
      (((c : ℕ+) : ℕ) * ((k : ℕ+) : ℕ)) with hlt | hge
  · exact hx₀ _ (by omega) (nsmul_sub_eq_zero_of_eq (le_of_lt hlt) hk')
  · have hne' : ((d : ℕ+) : ℕ) * ((k : ℕ+) : ℕ) ≠ ((c : ℕ+) : ℕ) * ((k : ℕ+) : ℕ) := by
      intro h
      exact hne (PNat.coe_injective (Nat.eq_of_mul_eq_mul_right hkpos h)).symm
    have hgt : ((c : ℕ+) : ℕ) * ((k : ℕ+) : ℕ) < ((d : ℕ+) : ℕ) * ((k : ℕ+) : ℕ) := by omega
    exact hx₀ _ (by omega) (nsmul_sub_eq_zero_of_eq (le_of_lt hgt) hk'.symm)

/-- ★★**第 2 条**(無限位数の単元)—— `x₀` が `Φ^birat` に入っていればよい。 -/
theorem birat_frobeniusCompact_cond2_of_mem (A : BiratCat P G)
    (x₀ : Φ.val (P.toElem.obj (biratDown P G A)).base)
    (hx₀mem : toGp _ x₀ ∈ phiBiratAt P G A)
    (hx₀ : ∀ n : ℕ, 0 < n →
      n • toGp (Φ.val (P.toElem.obj (biratDown P G A)).base) x₀ ≠ 0) :
    ∃ u : End A, u ∈ OTimes (biratPre P G) A ∧
      ∀ k : ℕ+, (u ^ ((k : ℕ+) : ℕ) : End A) ≠ 1 := by
  obtain ⟨u, hu, hux⟩ := hx₀mem
  refine ⟨u, hu, fun k hk => ?_⟩
  have h := congrArg (fun t : End A => biratDivGp ((t : A ⟶ A))) hk
  rw [biratDivGp_pow hu, hux] at h
  refine hx₀ ((k : ℕ+) : ℕ) k.2 ?_
  rw [h]
  exact biratDivGp_id A

/-- ★★★★★**`𝒞^birat` の Frobenius-compact 対象(逸脱 (B) を外した版)**。 -/
theorem birat_isFrobeniusCompact_of_primary
    (hfn : ∀ X : BiratCat P G, IsFrobeniusNormalized (biratPre P G) X)
    (hndOn : MonoidOn.IsNonDilatingOn Φ)
    (A : BiratCat P G)
    (hprim : ∀ x : Φ.val (P.toElem.obj (biratDown P G A)).base, IsPrimaryElt x →
      toGp _ x ∈ phiBiratAt P G A)
    (x₀ : Φ.val (P.toElem.obj (biratDown P G A)).base)
    (hx₀mem : toGp _ x₀ ∈ phiBiratAt P G A) (hx₀ : x₀ ≠ 0) :
    IsFrobeniusCompact (biratPre P G) A :=
  haveI hne := toGp_nsmul_ne_zero_of_ne_zero (P.divisorial _).1.1 (P.divisorial _).2 hx₀
  ⟨birat_frobeniusCompact_cond1 hfn A,
   birat_frobeniusCompact_cond2_of_mem P G A x₀ hx₀mem hne,
   fun θ c d hyp =>
     frobeniusCompact_cond3_of_eq (biratPre P G) θ c d
       (birat_cd_eq_of_primary P G A θ c d (hndOn _ (biratBase θ.inv)) hprim x₀ hx₀mem hne hyp)
       hyp⟩

/-- ★★★**`(𝒞^un-tr)^birat` の Frobenius-compact 対象(逸脱 (B) を外した版)**。 -/
theorem unTrBiratCompact_of_primary (G' : Frobenioid P)
    (hfn : ∀ X : BiratCat (unTrPre P Fc) (unTr_frobenioid P Fc G'),
      IsFrobeniusNormalized (unTrBiratPre P Fc G') X)
    (hndOn : MonoidOn.IsNonDilatingOn Φ)
    (A : BiratCat (unTrPre P Fc) (unTr_frobenioid P Fc G'))
    (hprim : ∀ x : Φ.val ((unTrPre P Fc).toElem.obj
        (biratDown (unTrPre P Fc) (unTr_frobenioid P Fc G') A)).base, IsPrimaryElt x →
      toGp _ x ∈ phiBiratAt (unTrPre P Fc) (unTr_frobenioid P Fc G') A)
    (x₀ : Φ.val ((unTrPre P Fc).toElem.obj
      (biratDown (unTrPre P Fc) (unTr_frobenioid P Fc G') A)).base)
    (hx₀mem : toGp _ x₀ ∈ phiBiratAt (unTrPre P Fc) (unTr_frobenioid P Fc G') A)
    (hx₀ : x₀ ≠ 0) :
    ∃ X : BiratCat (unTrPre P Fc) (unTr_frobenioid P Fc G'),
      IsFrobeniusCompact (unTrBiratPre P Fc G') X :=
  ⟨A, birat_isFrobeniusCompact_of_primary (unTrPre P Fc) (unTr_frobenioid P Fc G')
    hfn hndOn A hprim x₀ hx₀mem hx₀⟩

/-- ★★★locator —— `Theorem 6.2, (iii)` の Frobenius-compact 対象(逸脱 (B) なし)。 -/
def birat_isFrobeniusCompact_of_primary.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 111,
    item := "Theorem 6.2, (iii) — 𝒞^birat の Frobenius-compact 対象(台の条件だけで足りる)",
    sectionId := "frdi-thm-6-2" }

end NoDeviation

/-! ## ★5. model Frobenioid が rationally standard 型(逸脱 (B) なし) -/

namespace ModelData

variable {M : ModelData.{v, u, w} D}

/-- ★★★★★★★**[FrdI] Theorem 6.2, (iii)**(**逸脱 (B) なし**の版)——
model Frobenioid が **rationally standard 型**。

★★`unTrBiratCompact` の側を `birat_isFrobeniusCompact_of_primary` に差し替えたので、
必要なのは
* **準素元が `Φ^birat` に入ること**(原文の「台に入る」がちょうどこれ)
* 0 でない元が 1 つ `Φ^birat` に入ること(group-like でないこと)

だけであり、★**`Div : 𝒪^▷ ↠ Φ`(逸脱 (B))は要らない**。

原文 (FrdI p.111):
> of an element of B(L), then C is of rationally standard type. -/
theorem model_isOfRationallyStandardType_of_primary (h : Hyp M)
    (ι : ∀ Y : D, Prime (M.phi.val Y) → Pf (M.phi.val Y) → NNReal)
    (R : RatFnData (modelPre h) (model_frobenioid h))
    (hsp : ∀ (A : Obj M) (p : Prime (M.phi.val ((modelPre h).toElem.obj A).base)),
      ∃ (a b : M.phi.val ((modelPre h).toElem.obj A).base)
        (y : R.bmon.val ((modelPre h).toElem.obj A).base),
        (toGp _ a - toGp _ b = R.divB _ y ∨ toGp _ a - toGp _ b = -(R.divB _ y)) ∧
        p ∈ SuppElt (ι _) a ∧ p ∉ SuppElt (ι _) b)
    (hfsmff : IsOfFSMFFType D) (hnd : M.phi.IsNonDilatingOn)
    (hgl : IsOfGroupLikeType (modelPre h) →
      ∃ A, IsFrobeniusCompact (istrPre (modelPre h) (model_frobenioid h).core) A)
    (A : BiratCat (unTrPre (modelPre h) (model_frobenioid h).core)
      (unTr_frobenioid (modelPre h) (model_frobenioid h).core (model_frobenioid h)))
    (hprim : ∀ x : M.phi.val ((unTrPre (modelPre h) (model_frobenioid h).core).toElem.obj
        (biratDown _ _ A)).base, IsPrimaryElt x →
      toGp _ x ∈ phiBiratAt (unTrPre (modelPre h) (model_frobenioid h).core)
        (unTr_frobenioid (modelPre h) (model_frobenioid h).core (model_frobenioid h)) A)
    (x₀ : M.phi.val ((unTrPre (modelPre h) (model_frobenioid h).core).toElem.obj
      (biratDown _ _ A)).base)
    (hx₀mem : toGp _ x₀ ∈ phiBiratAt (unTrPre (modelPre h) (model_frobenioid h).core)
      (unTr_frobenioid (modelPre h) (model_frobenioid h).core (model_frobenioid h)) A)
    (hx₀ : x₀ ≠ 0) :
    haveI := h.connectedD
    IsOfRationallyStandardType (modelPre h) (model_frobenioid h) ι :=
  haveI := h.connectedD
  (model_ratStandardType_iff h ι).mpr
    ⟨isOfRationalType_of_divB R ι hsp,
     model_isOfStandardType h (model_frobenioid h).core hfsmff hnd hgl,
     unTrBiratCompact_of_primary (modelPre h) (model_frobenioid h).core (model_frobenioid h)
       (fun X => (unTr_isOfModelType (model_frobenioid h).core (model_frobenioid h)).2 X)
       hnd A hprim x₀ hx₀mem hx₀⟩

/-- ★★★★★★locator —— `Theorem 6.2, (iii)` の rationally standard 型(逸脱なし)。 -/
def model_isOfRationallyStandardType_of_primary.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 111,
    item := "Theorem 6.2, (iii) — model Frobenioid は rationally standard 型(逸脱 (B) なし)",
    sectionId := "frdi-thm-6-2" }

end ModelData

end ABC3.Found.FrdI
