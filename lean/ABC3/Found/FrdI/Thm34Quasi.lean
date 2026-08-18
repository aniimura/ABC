/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Thm34Pre
import ABC3.Found.FrdI.Prop33UnTr
import Mathlib.Data.PNat.Factors

/-!
# [FrdI] Theorem 3.4, (ii) を quasi-isotropic 型へ

★現行の `thm_3_4_ii` は isotropic 型を仮定しているが、
原典は quasi-isotropic 型で十分だと言う。その帰着の材料を置く。
-/

universe v u w u2 v2

namespace ABC3.Found.FrdI

open CategoryTheory

/-! ## ★★★`Theorem 3.4, (ii)` を quasi-isotropic 型へ一般化する準備

原文 (FrdI p.62):
> (ii) Suppose that C1, C2 are of quasi-isotropic type, and that D1, D2 are of

★★現行の `thm_3_4_ii` は **isotropic 型**を仮定しているが、
原典は **quasi-isotropic 型**で十分だと言う。
★帰着の筋は `𝒞^istr` へ落とすことで、
\`psiIstr_isEquivalence\`(quasi-isotropic を仮定として受ける)が既にある。

★★**欠けていた環**は「`𝒞` での pre-step 性 ⇔ `𝒞^istr` での pre-step 性」であり、
★在庫の `isotropification_degFr`(次数の保存)と
`isotropification_baseIso_iff`(底同型の両向き保存)を組むだけで出る。 -/

/-- ★★★★**isotropification は pre-step 性を両向きに保つ**。

★`ii-quasi`(quasi-isotropic への一般化)の欠けていた環。 -/
theorem isotropification_isPreStep_iff {Dd : Type u} [Category.{v} Dd] {Cc : Type u2}
    [Category.{v2} Cc] {Φ₀ : MonoidOn.{v, u, w} Dd} (P : PreFrobenioid Cc Φ₀)
    (F : FrobenioidCore P) {A B : Cc} (f : A ⟶ B) :
    IsPreStep P f ↔ IsPreStep (istrPre P F) ((isotropification P F).map f) := by
  constructor
  · rintro ⟨hl, hb⟩
    refine ⟨?_, ?_⟩
    · show P.degFr (istrMap P F f) = 1
      rw [isotropification_degFr]
      exact hl
    · exact (isotropification_baseIso_iff P F f).mpr hb
  · rintro ⟨hl, hb⟩
    refine ⟨?_, ?_⟩
    · show P.degFr f = 1
      rw [← isotropification_degFr P F f]
      exact hl
    · exact (isotropification_baseIso_iff P F f).mp hb

def isotropification_isPreStep_iff.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 62,
    item := "Theorem 3.4, (ii) — isotropification は pre-step を保つ",
    sectionId := "frdi-thm-3-4" }

/-! ## ★同型四角形で pre-step 性は不変

★★`isotropificationCommute`(`Thm34.lean:774`)の自然同型を渡るのに要る。
★`Ψ^istr ∘ istr₁` と `istr₂ ∘ Ψ` は**同型だが同一ではない**ので、
射の性質を四角形を越えて運ぶ補題が必要である。 -/

/-- ★★★★**同型の四角形を越えて pre-step 性は保たれる**。

★次数は同型が 1 だから、底同型は同型の合成から出る。 -/
theorem isPreStep_congr_iso {Dd : Type u} [Category.{v} Dd] {Cc : Type u2}
    [Category.{v2} Cc] {Φ₀ : MonoidOn.{v, u, w} Dd} (P : PreFrobenioid Cc Φ₀)
    {A B A' B' : Cc} (f : A ⟶ B) (f' : A' ⟶ B')
    (α : A ≅ A') (β : B ≅ B') (h : f ≫ β.hom = α.hom ≫ f') :
    IsPreStep P f ↔ IsPreStep P f' := by
  have hda : P.degFr α.hom = 1 := degFr_of_isIso P α.hom
  have hdb : P.degFr β.hom = 1 := degFr_of_isIso P β.hom
  have hd : P.degFr f = P.degFr f' := by
    have := congrArg P.degFr h
    rwa [P.degFr_comp, P.degFr_comp, hda, hdb, mul_one, one_mul] at this
  have hb := congrArg P.Base h
  rw [P.Base_comp, P.Base_comp] at hb
  haveI : IsIso (P.Base α.hom) := ⟨P.Base α.inv, by
    rw [← P.Base_comp, α.hom_inv_id, P.Base_id], by
    rw [← P.Base_comp, α.inv_hom_id, P.Base_id]⟩
  haveI : IsIso (P.Base β.hom) := ⟨P.Base β.inv, by
    rw [← P.Base_comp, β.hom_inv_id, P.Base_id], by
    rw [← P.Base_comp, β.inv_hom_id, P.Base_id]⟩
  constructor
  · rintro ⟨hl, hbi⟩
    refine ⟨?_, ?_⟩
    · show P.degFr f' = 1
      rw [← hd]; exact hl
    · show IsIso (P.Base f')
      haveI : IsIso (P.Base f) := hbi
      have : P.Base f' = inv (P.Base α.hom) ≫ (P.Base f ≫ P.Base β.hom) := by
        rw [hb, ← Category.assoc, IsIso.inv_hom_id, Category.id_comp]
      rw [this]; infer_instance
  · rintro ⟨hl, hbi⟩
    refine ⟨?_, ?_⟩
    · show P.degFr f = 1
      rw [hd]; exact hl
    · show IsIso (P.Base f)
      haveI : IsIso (P.Base f') := hbi
      have : P.Base f = (P.Base α.hom ≫ P.Base f') ≫ inv (P.Base β.hom) := by
        rw [← hb, Category.assoc, IsIso.hom_inv_id, Category.comp_id]
      rw [this]; infer_instance

def isPreStep_congr_iso.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 62,
    item := "Theorem 3.4, (ii) — 同型四角形で pre-step 性は不変",
    sectionId := "frdi-thm-3-4" }

/-! ## ★★★★`Theorem 3.4, (ii)` の quasi-isotropic 版

原文 (FrdI p.62):
> (ii) Suppose that C1, C2 are of quasi-isotropic type, and that D1, D2 are of

★★**帰着で出す**: `𝒞^istr` へ落として既存の isotropic 版を使う。 -/

variable {D₁ : Type u} [Category.{v} D₁] {C₁ : Type u2} [Category.{v2} C₁]
  {Φ₁ : MonoidOn.{v, u, w} D₁} {P₁ : PreFrobenioid C₁ Φ₁}
  {D₂ : Type u} [Category.{v} D₂] {C₂ : Type u2} [Category.{v2} C₂]
  {Φ₂ : MonoidOn.{v, u, w} D₂} {P₂ : PreFrobenioid C₂ Φ₂}

/-- ★★★★★**[FrdI] Theorem 3.4, (ii)** の **quasi-isotropic 版**(pre-step の条)。

★手: `isotropification_isPreStep_iff` で `𝒞^istr` へ落とし、
`thm_3_4_ii` を `Ψ^istr` に適用し、
`isotropificationCommute` の自然同型を `isPreStep_congr_iso` で渡る。 -/
theorem thm_3_4_ii_quasi_preStep (Ψ : C₁ ⥤ C₂) [Ψ.IsEquivalence]
    (F₁ : FrobenioidCore P₁) (G₁ : Frobenioid P₁)
    (hq₁ : IsOfQuasiIsotropicType C₁ P₁) (hFSMFF₁ : IsOfFSMFFType D₁)
    (hFrobMono₁ : ∀ {X Y : Istr P₁} (ε : X ⟶ Y),
      IsFrobeniusType (istrPre P₁ F₁) ε → Mono ε)
    (hFrobFS₁ : ∀ {X Y : Istr P₁} (ε : X ⟶ Y),
      IsFrobeniusType (istrPre P₁ F₁) ε → IsFiberwiseSurjective ε)
    (F₂ : FrobenioidCore P₂) (G₂ : Frobenioid P₂)
    (hq₂ : IsOfQuasiIsotropicType C₂ P₂) (hFSMFF₂ : IsOfFSMFFType D₂)
    (hFrobMono₂ : ∀ {X Y : Istr P₂} (ε : X ⟶ Y),
      IsFrobeniusType (istrPre P₂ F₂) ε → Mono ε)
    (hFrobFS₂ : ∀ {X Y : Istr P₂} (ε : X ⟶ Y),
      IsFrobeniusType (istrPre P₂ F₂) ε → IsFiberwiseSurjective ε)
    {A B : C₁} (φ : A ⟶ B) :
    IsPreStep P₁ φ ↔ IsPreStep P₂ (Ψ.map φ) := by
  haveI := psiIstr_isEquivalence Ψ P₁ P₂ hq₁ hq₂
  have key : ∀ {X Y : Istr P₁} (f : X ⟶ Y),
      IsPreStep (istrPre P₁ F₁) f ↔
        IsPreStep (istrPre P₂ F₂) ((psiIstr Ψ P₁ P₂ hq₁ hq₂).map f) :=
    fun f => (thm_3_4_ii (istrPre P₁ F₁) (istrPre P₂ F₂)
      (psiIstr Ψ P₁ P₂ hq₁ hq₂).asEquivalence
      (istr_frobenioidCore P₁ F₁) (istr_frobenioid P₁ F₁ G₁)
      (istr_isotropic P₁ F₁) hFSMFF₁ hFrobMono₁ hFrobFS₁
      (istr_frobenioidCore P₂ F₂) (istr_frobenioid P₂ F₂ G₂)
      (istr_isotropic P₂ F₂) hFSMFF₂ hFrobMono₂ hFrobFS₂).1 f
  have hcom := isotropificationCommute Ψ P₁ P₂ F₁ F₂ hq₁ hq₂
  rw [isotropification_isPreStep_iff P₁ F₁ φ, isotropification_isPreStep_iff P₂ F₂ (Ψ.map φ)]
  refine (key ((isotropification P₁ F₁).map φ)).trans ?_
  exact isPreStep_congr_iso (istrPre P₂ F₂) _ _ (hcom.app A) (hcom.app B)
    (hcom.hom.naturality φ)

def thm_3_4_ii_quasi_preStep.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 62,
    item := "Theorem 3.4, (ii) — quasi-isotropic 型での pre-step の保存",
    sectionId := "frdi-thm-3-4" }

/-! ## ★★素数の全単射から `ℕ≥1 ≃* ℕ≥1` を作る

★★`Theorem 3.4, (iii)` の `Ψ_{N≥1}` の**拡張段**である。

原文 (FrdI p.62):
> such that Ψ maps morphisms of Frobenius degree d to morphisms of Frobenius

★`ℕ≥1` は素数の自由可換単系なので、素数の全単射がそのまま単系自己同型を与える。
★mathlib の `PNat.factorMultiset` / `PrimeMultiset` を使う。 -/

/-- ★★★★**素数の全単射 ⇒ `ℕ≥1` の単系自己同型**。 -/
def pnatMulEquivOfPrimeEquiv (σ : Nat.Primes ≃ Nat.Primes) : ℕ+ ≃* ℕ+ where
  toFun n := PrimeMultiset.prod (Multiset.map σ (PNat.factorMultiset n))
  invFun n := PrimeMultiset.prod (Multiset.map σ.symm (PNat.factorMultiset n))
  left_inv n := by
    show PrimeMultiset.prod (Multiset.map σ.symm (PNat.factorMultiset
      (PrimeMultiset.prod (Multiset.map σ (PNat.factorMultiset n))))) = n
    rw [PrimeMultiset.factorMultiset_prod, Multiset.map_map,
      show (σ.symm ∘ σ) = id from funext σ.symm_apply_apply, Multiset.map_id]
    exact PNat.prod_factorMultiset n
  right_inv n := by
    show PrimeMultiset.prod (Multiset.map σ (PNat.factorMultiset
      (PrimeMultiset.prod (Multiset.map σ.symm (PNat.factorMultiset n))))) = n
    rw [PrimeMultiset.factorMultiset_prod, Multiset.map_map,
      show (σ ∘ σ.symm) = id from funext σ.apply_symm_apply, Multiset.map_id]
    exact PNat.prod_factorMultiset n
  map_mul' m n := by
    show PrimeMultiset.prod (Multiset.map σ (PNat.factorMultiset (m * n)))
      = PrimeMultiset.prod (Multiset.map σ (PNat.factorMultiset m))
        * PrimeMultiset.prod (Multiset.map σ (PNat.factorMultiset n))
    rw [PNat.factorMultiset_mul]
    show PrimeMultiset.prod (Multiset.map σ
        (PNat.factorMultiset m + PNat.factorMultiset n)) = _
    have hma : Multiset.map σ (PNat.factorMultiset m + PNat.factorMultiset n)
        = Multiset.map σ (PNat.factorMultiset m) + Multiset.map σ (PNat.factorMultiset n) :=
      Multiset.map_add _ _ _
    rw [hma]
    exact PrimeMultiset.prod_add _ _

def pnatMulEquivOfPrimeEquiv.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 62,
    item := "Theorem 3.4, (iii) — Ψ_{N≥1} の拡張段",
    sectionId := "frdi-thm-3-4" }

/-! ## ★`psiN-prime-map` の部品

★★`primeFrobDivId_map_iff`(`iii-primefrob-full`)は
`IsPrimeFrobenius ∧ IsDivIdentity` を運ぶが、
`IsDivIdentity` は**自己射のみ**に定義されている。
★★★**測定の訂正(2026-08-18)**: 一度「Frobenius 型なら `IsDivIdentity` は自動」と
書こうとしたが、**偽である**。
★FrdI の `Div-equivalent` は `Φ(φ) = Φ(ψ)`(誘導される単項式射の一致)であって、
`metrically equivalent`(`Div(φ) = Div(ψ)`)とは**別条件**である
(`MorphismTypes.lean:148` の注意書きのとおり)。
★★したがって `IsDivIdentity` は `Φ.map (Base φ) = Φ.map (𝟙)` であり、
**`Φ` が non-dilating であること**から出る筋になる。
★★★これが `Theorem 3.4, (iii)` が **standard 型**を仮定する理由である
(`IsOfStandardType` の (e) が `Φ` の non-dilating)。 -/

/-- ★★**素数次数の Frobenius 型射は prime-Frobenius**。 -/
theorem isPrimeFrobenius_of_degFr_prime {Dd : Type u} [Category.{v} Dd] {Cc : Type u2}
    [Category.{v2} Cc] {Φ₀ : MonoidOn.{v, u, w} Dd} (P : PreFrobenioid Cc Φ₀)
    {A B : Cc} (φ : A ⟶ B) (h : IsFrobeniusType P φ) {p : ℕ+}
    (hd : P.degFr φ = p) (hp : Nat.Prime ((p : ℕ+) : ℕ)) :
    IsPrimeFrobenius P φ :=
  ⟨h, by rw [hd]; exact hp⟩

def isPrimeFrobenius_of_degFr_prime.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 62,
    item := "Theorem 3.4, (iii) — 素数次数の Frobenius 型射",
    sectionId := "frdi-thm-3-4" }

/-! ## ★★★`(p₁, p₂)-admissible` —— 原典の (iii) の骨格

原文 (FrdI p.64):
> (p1, p2)-admissible if Ψ maps every p1-Frobenius morphism with domain A1 to a

★★原典は素数の対応を次の 2 つで出す:
- **(F1)** 各素数 `p₁` に対し、ある素数 `p₂` と `(p₁,p₂)-admissible` な対象がある。
- **(F2)** 任意の射 `ζ₁ : A₁ ⟶ B₁` に対し、`A₁` が admissible ⇔ `B₁` が admissible。

★そして **`𝒞₁` の連結性**から全対象が admissible になる。 -/

variable {D₃ : Type u} [Category.{v} D₃] {C₃ : Type u2} [Category.{v2} C₃]
  {Φ₃ : MonoidOn.{v, u, w} D₃} {P₃ : PreFrobenioid C₃ Φ₃}
  {D₄ : Type u} [Category.{v} D₄] {C₄ : Type u2} [Category.{v2} C₄]
  {Φ₄ : MonoidOn.{v, u, w} D₄} {P₄ : PreFrobenioid C₄ Φ₄}

/-- ★★★★**[FrdI] Theorem 3.4, (iii)** の `(p₁, p₂)-admissible`。

`A₁` から出るすべての `p₁`-Frobenius 射を `Ψ` が `p₂`-Frobenius 射へ送ること。 -/
def IsAdmissibleObj (P₁ : PreFrobenioid C₃ Φ₃) (P₂ : PreFrobenioid C₄ Φ₄)
    (Ψ : C₃ ⥤ C₄) (p₁ p₂ : ℕ+) (A : C₃) : Prop :=
  ∀ {B : C₃} (φ : A ⟶ B), IsFrobeniusType P₁ φ → P₁.degFr φ = p₁ →
    IsFrobeniusType P₂ (Ψ.map φ) ∧ P₂.degFr (Ψ.map φ) = p₂

/-- ★★★**連結性で (F1)＋(F2) から全対象へ広がる**。

★原典が (F1)(F2) を置いた後に使う歩みを、射の連の形で書いた。
★実際の連結性(`P.connectedC`)は zigzag を与えるので、
この補題を zigzag に沿って繰り返すことになる。 -/
theorem isAdmissibleObj_of_hom (P₁ : PreFrobenioid C₃ Φ₃) (P₂ : PreFrobenioid C₄ Φ₄)
    (Ψ : C₃ ⥤ C₄) (p₁ p₂ : ℕ+)
    (hF2 : ∀ {X Y : C₃} (_ : X ⟶ Y),
      IsAdmissibleObj P₁ P₂ Ψ p₁ p₂ X ↔ IsAdmissibleObj P₁ P₂ Ψ p₁ p₂ Y)
    {A B : C₃} (ζ : A ⟶ B) (hA : IsAdmissibleObj P₁ P₂ Ψ p₁ p₂ A) :
    IsAdmissibleObj P₁ P₂ Ψ p₁ p₂ B :=
  (hF2 ζ).mp hA

def IsAdmissibleObj.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 64,
    item := "Theorem 3.4, (iii) — (p₁,p₂)-admissible 対象",
    sectionId := "frdi-thm-3-4" }

/-! ## ★★★base-identity ⇒ Div-identity

原文 (FrdI p.65):
> exists a base-identity [hence Div-identity] p1-Frobenius endomorphism φ

★★★**この括弧書きが鍵である**。一度 `IsDivIdentity` を
**Frobenius 型**から出そうとして詰まったが、原典は
**base-identity** から出している。
★`IsDivIdentity φ = Φ.map (Base φ) = Φ.map (Base 𝟙)` なので、
`Base φ = Base 𝟙` に `congrArg` を当てるだけで出る。
★実は**既に在庫にあった** —— `MorphismTypes.lean:530` の `isDivIdentity_of_isBaseIdentity`。 -/

/-! ## ★★★(F1) の group-like の場合

原文 (FrdI p.65):
> with domain A1 are precisely the irreducible base-isomorphisms with domain A1 [cf.

★★group-like 型では **prime-Frobenius 射 ＝ 既約な底同型**である。
★右辺は**圧倒的な条件だけ**で書けるので、`Ψ` が保つことが直ちに出る。 -/

/-- ★★★★**group-like 型では prime-Frobenius ⇔ 既約な底同型**。

★手: `prop_1_14_i` の三分法で、
- step の枝は group-like だから pre-step が同型になり矛盾
- pull-back の枝は `Base φ` が同型なので既約になりえず矛盾 -/
theorem isPrimeFrobenius_iff_irred_baseIso_of_groupLike
    {Dd : Type u} [Category.{v} Dd] {Cc : Type u2} [Category.{v2} Cc]
    {Φ₀ : MonoidOn.{v, u, w} Dd} (P : PreFrobenioid Cc Φ₀) (G : Frobenioid P)
    (hiso : ∀ X : Cc, IsIsotropic P X) (hgl : IsOfGroupLikeType P)
    {A B : Cc} (φ : A ⟶ B) :
    IsPrimeFrobenius P φ ↔ (IsIrreducibleMor φ ∧ IsBaseIsomorphism P φ) := by
  constructor
  · intro h
    exact ⟨prop_1_10_iv_mp P G.core (hiso A) φ h, h.1.2⟩
  · rintro ⟨hirr, hb⟩
    rcases (prop_1_14_i P G hiso φ).mp hirr with hpf | ⟨hst, _⟩ | ⟨_, hbirr⟩
    · exact hpf
    · exact absurd (isIso_of_preStep_of_isGroupLikeObj P G.core
        (fun X _ => hiso X) (hgl A) φ hst.1) hst.2
    · haveI : IsIso (P.Base φ) := hb
      exact absurd hbirr (by
        simpa using not_isIrreducibleMor_of_isIso (P.Base φ))

def isPrimeFrobenius_iff_irred_baseIso_of_groupLike.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 65,
    item := "Theorem 3.4, (iii) (F1) — group-like での prime-Frobenius の圧倒的特徴づけ",
    sectionId := "frdi-thm-3-4" }

/-! ## ★(F1) の (b) —— Frobenius-trivial 対象の base-identity な p-Frobenius 自己射

原文 (FrdI p.65):
> exists a base-identity [hence Div-identity] p1-Frobenius endomorphism φ

★★**定義から直ちに出る**。`IsFrobeniusTrivial A` は
「各次数の base-identity な Frobenius 型自己射の族」そのものだからである。 -/

/-- ★★★**Frobenius-trivial 対象には各素数次数の
base-identity な prime-Frobenius 自己射がある**。 -/
theorem exists_baseId_primeFrob_of_frobTrivial {Dd : Type u} [Category.{v} Dd]
    {Cc : Type u2} [Category.{v2} Cc] {Φ₀ : MonoidOn.{v, u, w} Dd}
    (P : PreFrobenioid Cc Φ₀) {A : Cc} (h : IsFrobeniusTrivial P A)
    {p : ℕ+} (hp : Nat.Prime ((p : ℕ+) : ℕ)) :
    ∃ φ : A ⟶ A, IsPrimeFrobenius P φ ∧ IsBaseIdentity P φ ∧ P.degFr φ = p := by
  obtain ⟨ζ, hdeg, hprop⟩ := h
  refine ⟨((ζ p : End A) : A ⟶ A), ⟨(hprop p).2, ?_⟩, (hprop p).1, hdeg p⟩
  rw [hdeg p]
  exact hp

def exists_baseId_primeFrob_of_frobTrivial.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 65,
    item := "Theorem 3.4, (iii) (F1) — Frobenius-trivial 対象の p-Frobenius 自己射",
    sectionId := "frdi-thm-3-4" }

/-! ## ★★「1 本から全部へ」—— `frobDegUniq` が効く

★★(F1) の admissible 性は「**すべての** p₁-Frobenius 射」についての主張だが、
原典は base-identity な自己射 **1 本**しか作らない。
★隔たりを埋めるのが `Definition 1.3, (ii)` の **`frobDegUniq`**
(同じ対象から出る同次数の Frobenius 型射は同型を除いて一意)である。 -/

/-- ★★★★**同じ対象から出る同次数の Frobenius 型射は、像の次数も等しい**。

★これで「1 本で判定すれば全部に及ぶ」となる。 -/
theorem degFr_map_eq_of_sameDeg
    (F₁ : FrobenioidCore P₁) (Ψ : C₁ ⥤ C₂)
    {A B E : C₁} (φ : A ⟶ B) (ψ : A ⟶ E)
    (hφ : IsFrobeniusType P₁ φ) (hψ : IsFrobeniusType P₁ ψ)
    (hd : P₁.degFr φ = P₁.degFr ψ) :
    P₂.degFr (Ψ.map φ) = P₂.degFr (Ψ.map ψ) := by
  obtain ⟨β, hβiso, hβ⟩ := F₁.frobDegUniq A B E φ ψ hφ hψ hd
  haveI := hβiso
  haveI : IsIso (Ψ.map β) := Ψ.map_isIso β
  have h : P₂.degFr (Ψ.map (φ ≫ β)) = P₂.degFr (Ψ.map ψ) := by rw [hβ]
  rw [Ψ.map_comp, P₂.degFr_comp, degFr_of_isIso P₂ (Ψ.map β), one_mul] at h
  exact h

def degFr_map_eq_of_sameDeg.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 65,
    item := "Theorem 3.4, (iii) (F1) — 同次数の Frobenius 型射は像の次数も等しい",
    sectionId := "frdi-thm-3-4" }

/-! ## ★`Proposition 1.8, (iii)` の第 3 主張は既に在庫にあった

原文 (FrdI p.30):
> base-isomorphic objects, then A is group-like if and only if B is.

★★**在庫検索の失敗を本セッションで 5 回目として記録する**。
実装しようとしたところ Lean が重複を報告した ——
`Prop110.lean:2733` の `isGroupLikeObj_of_baseIso`
(`α : A' ⟶ A` と `IsBaseIsomorphism P α` を取る形)。

★★★これで **(F1) の (a) の部品は全部揃った**:
- `baseSurj`(= `Definition 1.3, (i)(a)`)が Frobenius-trivial 対象を与える
- `isGroupLikeObj_of_baseIso` が非 group-like 性を運ぶ
- `preStepSpan` が底同型を pre-step の span に直す
-/

/-- ★★★★**group-like 性は底の同型を越えて両向きに移る**。

★在庫の `isGroupLikeObj_of_baseIso` は `𝒞` の射を取るので**片向き**だが、
`IsGroupLikeObj P A := IsGroupLike (Φ.val (Base A))` は**底にしか依らない**ので、
底の同型から両向きに出る。 -/
theorem isGroupLikeObj_of_baseIsoD {Dd : Type u} [Category.{v} Dd] {Cc : Type u2}
    [Category.{v2} Cc] {Φ₀ : MonoidOn.{v, u, w} Dd} (P : PreFrobenioid Cc Φ₀)
    {A B : Cc} (α : (P.toElem.obj A).base ≅ (P.toElem.obj B).base)
    (h : IsGroupLikeObj P A) : IsGroupLikeObj P B := by
  rw [isGroupLike_iff] at h ⊢
  intro b
  have hround : Φ₀.map α.inv (Φ₀.map α.hom b) = b := by
    have hc := Φ₀.map_comp α.hom α.inv b
    rw [α.inv_hom_id] at hc
    rw [hc]
    exact MonoidOn.map_id Φ₀ _ b
  rw [← hround]
  exact (h (Φ₀.map α.hom b)).map (Φ₀.map α.inv)

def isGroupLikeObj_of_baseIsoD.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 30,
    item := "Proposition 1.8, (iii) — 底の同型を越えた group-like 性",
    sectionId := "frdi-prop-1-8" }

end ABC3.Found.FrdI
