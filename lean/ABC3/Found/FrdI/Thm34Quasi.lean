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
  show IsGroupLike (Φ₀.val (P.toElem.obj B).base)
  have h : IsGroupLike (Φ₀.val (P.toElem.obj A).base) := h
  rw [isGroupLike_iff] at h ⊢
  intro b
  have hround : Φ₀.map α.inv (Φ₀.map α.hom b) = b := by
    have hc := Φ₀.map_comp α.hom α.inv b
    rw [α.inv_hom_id] at hc
    rw [← hc]
    exact MonoidOn.map_id Φ₀ _ b
  rw [← hround]
  exact (h (Φ₀.map α.hom b)).map (Φ₀.map α.inv)

def isGroupLikeObj_of_baseIsoD.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 30,
    item := "Proposition 1.8, (iii) — 底の同型を越えた group-like 性",
    sectionId := "frdi-prop-1-8-iii" }

/-! ## ★★★★(F1) の非 group-like の場合を束ねる

原文 (FrdI p.64):
> (p1, p2)-admissible object of C1.

★部品は全部揃っているので、ここで束ねる。 -/

/-- ★★★★**Frobenius-trivial かつ非 group-like な対象は admissible**。

★(F1) の非 group-like の場合そのもの。 -/
theorem exists_admissible_of_frobTrivial (e : C₁ ≌ C₂)
    (F₁ : FrobenioidCore P₁) (G₁ : Frobenioid P₁)
    (hiso₁ : ∀ X : C₁, IsIsotropic P₁ X) (hnd₁ : MonoidOn.IsNonDilatingOn Φ₁)
    (F₂ : FrobenioidCore P₂) (G₂ : Frobenioid P₂)
    (hiso₂ : ∀ X : C₂, IsIsotropic P₂ X) (hnd₂ : MonoidOn.IsNonDilatingOn Φ₂)
    (hps : ∀ {X Y : C₁} (f : X ⟶ Y), IsPreStep P₁ f ↔ IsPreStep P₂ (e.functor.map f))
    (hps' : ∀ {X Y : C₂} (f : X ⟶ Y), IsPreStep P₂ f ↔ IsPreStep P₁ (e.inverse.map f))
    {A : C₁} (hft : IsFrobeniusTrivial P₁ A)
    (hA : ¬ IsGroupLikeObj P₁ A) (hA₂ : ¬ IsGroupLikeObj P₂ (e.functor.obj A))
    {p₁ : ℕ+} (hp₁ : Nat.Prime ((p₁ : ℕ+) : ℕ)) :
    ∃ p₂ : ℕ+, IsAdmissibleObj P₁ P₂ e.functor p₁ p₂ A := by
  obtain ⟨φ₁, hpf₁, hbid₁, hdeg₁⟩ := exists_baseId_primeFrob_of_frobTrivial P₁ hft hp₁
  have hirr₁ : IsIrreducibleMor φ₁ := prop_1_10_iv_mp P₁ F₁ (hiso₁ A) φ₁ hpf₁
  have hmap := primeFrob_baseId_map P₁ P₂ e F₁ G₁ hiso₁ hnd₁ F₂ G₂ hiso₂ hnd₂
    hps hps' hA hA₂ φ₁ hirr₁ hpf₁ hbid₁
  refine ⟨P₂.degFr (e.functor.map φ₁), ?_⟩
  intro B ψ hψ hdψ
  have hsame : P₁.degFr φ₁ = P₁.degFr ψ := by rw [hdeg₁, hdψ]
  obtain ⟨β, hβiso, hβ⟩ := F₁.frobDegUniq A _ B φ₁ ψ hpf₁.1 hψ hsame
  haveI := hβiso
  refine ⟨?_, ?_⟩
  · rw [← hβ, e.functor.map_comp]
    exact IsFrobeniusType.comp P₂ F₂ hmap.1
      (isFrobeniusType_of_isIso P₂ (e.functor.map β))
  · exact (degFr_map_eq_of_sameDeg F₁ e.functor φ₁ ψ hpf₁.1 hψ hsame).symm

def exists_admissible_of_frobTrivial.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 65,
    item := "Theorem 3.4, (iii) (F1) — 非 group-like の場合",
    sectionId := "frdi-thm-3-4" }

/-- ★★★★**(F1) の group-like の場合** —— 任意の対象が admissible。

★`psiN-F1-gl`(prime-Frobenius ＝ 既約な底同型)と
`degFr_map_eq_of_sameDeg`(1 本から全部へ)を組む。 -/
theorem exists_admissible_of_groupLike (e : C₁ ≌ C₂)
    (F₁ : FrobenioidCore P₁) (G₁ : Frobenioid P₁)
    (hiso₁ : ∀ X : C₁, IsIsotropic P₁ X) (hgl₁ : IsOfGroupLikeType P₁)
    (G₂ : Frobenioid P₂)
    (hiso₂ : ∀ X : C₂, IsIsotropic P₂ X) (hgl₂ : IsOfGroupLikeType P₂)
    (hbiso : ∀ {X Y : C₁} (f : X ⟶ Y),
      IsBaseIsomorphism P₁ f → IsBaseIsomorphism P₂ (e.functor.map f))
    (A : C₁) {p₁ : ℕ+} (hp₁ : Nat.Prime ((p₁ : ℕ+) : ℕ)) :
    ∃ p₂ : ℕ+, IsAdmissibleObj P₁ P₂ e.functor p₁ p₂ A := by
  obtain ⟨B₀, ψ₀, hft₀, hd₀⟩ := F₁.frobDegSurj A p₁
  have hpf₀ : IsPrimeFrobenius P₁ ψ₀ := ⟨hft₀, by rw [hd₀]; exact hp₁⟩
  refine ⟨P₂.degFr (e.functor.map ψ₀), ?_⟩
  intro B ψ hψ hdψ
  have hsame : P₁.degFr ψ₀ = P₁.degFr ψ := by rw [hd₀, hdψ]
  have hpf : IsPrimeFrobenius P₁ ψ := ⟨hψ, by rw [hdψ]; exact hp₁⟩
  have hirr : IsIrreducibleMor ψ :=
    ((isPrimeFrobenius_iff_irred_baseIso_of_groupLike P₁ G₁ hiso₁ hgl₁ ψ).mp hpf).1
  refine ⟨?_, (degFr_map_eq_of_sameDeg F₁ e.functor ψ₀ ψ hft₀ hψ hsame).symm⟩
  exact ((isPrimeFrobenius_iff_irred_baseIso_of_groupLike P₂ G₂ hiso₂ hgl₂
    (e.functor.map ψ)).mpr
      ⟨isIrreducibleMor_map e.functor hirr, hbiso ψ hψ.2⟩).1

def exists_admissible_of_groupLike.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 64,
    item := "Theorem 3.4, (iii) (F1) — group-like の場合",
    sectionId := "frdi-thm-3-4" }

/-! ## ★★★★★(F1) 全体

原文 (FrdI p.64):
> (p1, p2)-admissible object of C1. -/

/-- ★★★★★**[FrdI] Theorem 3.4, (iii) の (F1)** ——
**各素数 `p₁` に対し、ある素数 `p₂` と `(p₁,p₂)-admissible` な対象がある**。

★group-like か否かで場分けする。非 group-like 側では
`baseSurj` で Frobenius-trivial な対象を作り、
`isGroupLikeObj_of_baseIsoD` で非 group-like 性を運ぶ。 -/
theorem exists_admissible_F1 (e : C₁ ≌ C₂)
    (F₁ : FrobenioidCore P₁) (G₁ : Frobenioid P₁)
    (hiso₁ : ∀ X : C₁, IsIsotropic P₁ X) (hnd₁ : MonoidOn.IsNonDilatingOn Φ₁)
    (F₂ : FrobenioidCore P₂) (G₂ : Frobenioid P₂)
    (hiso₂ : ∀ X : C₂, IsIsotropic P₂ X) (hnd₂ : MonoidOn.IsNonDilatingOn Φ₂)
    (hps : ∀ {X Y : C₁} (f : X ⟶ Y), IsPreStep P₁ f ↔ IsPreStep P₂ (e.functor.map f))
    (hps' : ∀ {X Y : C₂} (f : X ⟶ Y), IsPreStep P₂ f ↔ IsPreStep P₁ (e.inverse.map f))
    (hbiso : ∀ {X Y : C₁} (f : X ⟶ Y),
      IsBaseIsomorphism P₁ f → IsBaseIsomorphism P₂ (e.functor.map f))
    (hglmap : ∀ X : C₁, IsGroupLikeObj P₁ X ↔ IsGroupLikeObj P₂ (e.functor.obj X))
    (hglty : IsOfGroupLikeType P₁ → IsOfGroupLikeType P₂)
    (A₀ : C₁) {p₁ : ℕ+} (hp₁ : Nat.Prime ((p₁ : ℕ+) : ℕ)) :
    ∃ (p₂ : ℕ+) (A : C₁), IsAdmissibleObj P₁ P₂ e.functor p₁ p₂ A := by
  by_cases hgl : ∀ X : C₁, IsGroupLikeObj P₁ X
  · obtain ⟨p₂, hp₂⟩ := exists_admissible_of_groupLike e F₁ G₁ hiso₁ hgl G₂ hiso₂
      (hglty hgl) hbiso A₀ hp₁
    exact ⟨p₂, A₀, hp₂⟩
  · push Not at hgl
    obtain ⟨X₀, hX₀⟩ := hgl
    obtain ⟨A, hft, ⟨α⟩⟩ := F₁.baseSurj (P₁.toElem.obj X₀).base
    have hA : ¬ IsGroupLikeObj P₁ A := fun h =>
      hX₀ (isGroupLikeObj_of_baseIsoD P₁ α h)
    obtain ⟨p₂, hp₂⟩ := exists_admissible_of_frobTrivial e F₁ G₁ hiso₁ hnd₁
      F₂ G₂ hiso₂ hnd₂ hps hps' hft hA (fun h => hA ((hglmap A).mpr h)) hp₁
    exact ⟨p₂, A, hp₂⟩

def exists_admissible_F1.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 64,
    item := "Theorem 3.4, (iii) — (F1)",
    sectionId := "frdi-thm-3-4" }

/-! ## ★★★(F2) の次数移送

原文 (FrdI p.65):
> applying Proposition 1.14, (iv), to commutative diagrams such as the one given in

★★原典は `Proposition 1.11, (iii)` で四角形
`ψ₁ ∘ ζ₁ = ζ₁ ∘ φ₁` を作り、像でも四角形が保たれることから
`Proposition 1.14, (iv)` で次数を移す。
★在庫の `prop_1_14_iv_degFr` がそのまま使える。 -/

/-- ★★★★**四角形を越えて次数が移る** —— (F2) の中核。

`ζ ≫ φ = ψ ≫ ζ` なら、像でも `degFr (Ψφ) = degFr (Ψψ)`。 -/
theorem degFr_map_eq_of_square (Ψ : C₁ ⥤ C₂)
    {A B : C₁} (ζ : B ⟶ A) (φ : A ⟶ A) (ψ : B ⟶ B)
    (hsq : ζ ≫ φ = ψ ≫ ζ) :
    P₂.degFr (Ψ.map φ) = P₂.degFr (Ψ.map ψ) := by
  have hsq' : Ψ.map ζ ≫ Ψ.map φ = Ψ.map ψ ≫ Ψ.map ζ := by
    rw [← Ψ.map_comp, ← Ψ.map_comp, hsq]
  exact prop_1_14_iv_degFr P₂ (Ψ.map ζ) (Ψ.map φ) (Ψ.map ψ) (Ψ.map ζ) hsq' rfl

def degFr_map_eq_of_square.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 65,
    item := "Theorem 3.4, (iii) (F2) — 四角形を越えた次数の移送",
    sectionId := "frdi-thm-3-4" }

/-- ★★★★**(F2) の伝搬本体** —— 四角形があれば admissible 性が移る。

★`A` が admissible で、`ζ ≫ φ = ψ ≫ ζ` となる同次数の Frobenius 型自己射があれば、
`B` も admissible である。 -/
theorem admissible_of_square (Ψ : C₁ ⥤ C₂) (F₁ : FrobenioidCore P₁)
    {p₁ p₂ : ℕ+} {A B : C₁} (ζ : B ⟶ A) (φ : A ⟶ A) (ψ : B ⟶ B)
    (hsq : ζ ≫ φ = ψ ≫ ζ)
    (hφ : IsFrobeniusType P₁ φ) (hdφ : P₁.degFr φ = p₁)
    (hψ : IsFrobeniusType P₁ ψ) (hdψ : P₁.degFr ψ = p₁)
    (hFT : ∀ {X Y : C₁} (f : X ⟶ Y),
      IsFrobeniusType P₁ f → IsFrobeniusType P₂ (Ψ.map f))
    (hA : IsAdmissibleObj P₁ P₂ Ψ p₁ p₂ A) :
    IsAdmissibleObj P₁ P₂ Ψ p₁ p₂ B := by
  intro E θ hθ hdθ
  refine ⟨hFT θ hθ, ?_⟩
  have h1 : P₂.degFr (Ψ.map ψ) = P₂.degFr (Ψ.map θ) :=
    degFr_map_eq_of_sameDeg F₁ Ψ ψ θ hψ hθ (by rw [hdψ, hdθ])
  have h2 : P₂.degFr (Ψ.map φ) = P₂.degFr (Ψ.map ψ) :=
    degFr_map_eq_of_square Ψ ζ φ ψ hsq
  have h3 : P₂.degFr (Ψ.map φ) = p₂ := (hA φ hφ hdφ).2
  rw [← h1, ← h2, h3]

def admissible_of_square.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 65,
    item := "Theorem 3.4, (iii) (F2) — 四角形に沿う admissible 性の伝搬",
    sectionId := "frdi-thm-3-4" }

/-! ## ★★★(F2) の四角形を作る

原文 (FrdI p.65):
> that for every p1 ∈Primes, there exist base-identity p1-Frobenius endomorphisms

★★`Proposition 1.11, (iii)`(在庫、`prop_1_11_iii`)を
**base-identity** な自己射に当てると、`βD = 𝟙` で四角形が作れる。 -/

/-- ★★★★**pull-back に沿って base-identity 自己射が四角形を作る**。

★次数も一致する(`prop_1_14_iv_degFr`)。 -/
theorem exists_square_of_pullBack
    {A B : C₁} (ζ : B ⟶ A) (hpb : IsPullBack P₁ ζ)
    (α : End A) (hbid : IsBaseIdentity P₁ ((α : A ⟶ A))) :
    ∃ β : End B, IsBaseIdentity P₁ ((β : B ⟶ B)) ∧
      ζ ≫ ((α : A ⟶ A)) = ((β : B ⟶ B)) ≫ ζ ∧
      P₁.degFr ((β : B ⟶ B)) = P₁.degFr ((α : A ⟶ A)) := by
  have hbα : P₁.Base ((α : A ⟶ A)) = 𝟙 _ := by
    have h : P₁.Base ((α : A ⟶ A)) = P₁.Base (𝟙 A) := hbid
    rwa [P₁.Base_id] at h
  have hsq : P₁.Base ζ ≫ P₁.Base ((α : A ⟶ A))
      = ((1 : End ((P₁.toElem.obj B).base)) : _ ⟶ _) ≫ P₁.Base ζ := by
    rw [hbα, Category.comp_id]
    show P₁.Base ζ = 𝟙 _ ≫ P₁.Base ζ
    rw [Category.id_comp]
  obtain ⟨β, ⟨hbβ, hsqβ⟩, -⟩ := prop_1_11_iii P₁ ζ hpb α 1 hsq
  refine ⟨β, ?_, hsqβ, ?_⟩
  · show P₁.Base ((β : B ⟶ B)) = P₁.Base (𝟙 B)
    rw [P₁.Base_id, hbβ]
    rfl
  · exact (prop_1_14_iv_degFr P₁ ζ ((α : A ⟶ A)) ((β : B ⟶ B)) ζ hsqβ rfl).symm

def exists_square_of_pullBack.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 65,
    item := "Theorem 3.4, (iii) (F2) — pull-back に沿う四角形",
    sectionId := "frdi-thm-3-4" }

/-- ★★★★**(F2) の pre-step の場合**。

★`prop_1_10_i_exists_preStep`(在庫)が四角形 `ζ ≫ β = α ≫ φ'` を与え、
両辺の pre-step は次数 1 なので `prop_1_14_iv_degFr` がそのまま効く。 -/
theorem admissible_of_preStep (Ψ : C₁ ⥤ C₂) (F₁ : FrobenioidCore P₁)
    (hps : ∀ {X Y : C₁} (f : X ⟶ Y), IsPreStep P₁ f → IsPreStep P₂ (Ψ.map f))
    (hFT : ∀ {X Y : C₁} (f : X ⟶ Y),
      IsFrobeniusType P₁ f → IsFrobeniusType P₂ (Ψ.map f))
    {p₁ p₂ : ℕ+} {A B : C₁} (ζ : A ⟶ B) (hζ : IsPreStep P₁ ζ)
    (hA : IsAdmissibleObj P₁ P₂ Ψ p₁ p₂ A)
    {A' : C₁} (α : A ⟶ A') (hα : IsFrobeniusType P₁ α) (hdα : P₁.degFr α = p₁) :
    IsAdmissibleObj P₁ P₂ Ψ p₁ p₂ B := by
  obtain ⟨B', β, φ', hβ, hdβ, hφ', hsq⟩ :=
    prop_1_10_i_exists_preStep P₁ F₁ ζ hζ α hα
  intro E θ hθ hdθ
  refine ⟨hFT θ hθ, ?_⟩
  have h1 : P₂.degFr (Ψ.map β) = P₂.degFr (Ψ.map θ) :=
    degFr_map_eq_of_sameDeg F₁ Ψ β θ hβ hθ (by rw [hdβ, hdα, hdθ])
  have h2 : P₂.degFr (Ψ.map β) = P₂.degFr (Ψ.map α) := by
    have hsq' : Ψ.map ζ ≫ Ψ.map β = Ψ.map α ≫ Ψ.map φ' := by
      rw [← Ψ.map_comp, ← Ψ.map_comp, hsq]
    refine prop_1_14_iv_degFr P₂ (Ψ.map ζ) (Ψ.map β) (Ψ.map α) (Ψ.map φ') hsq' ?_
    rw [show P₂.degFr (Ψ.map ζ) = 1 from (hps ζ hζ).1,
      show P₂.degFr (Ψ.map φ') = 1 from (hps φ' hφ').1]
  rw [← h1, h2]
  exact (hA α hα hdα).2

def admissible_of_preStep.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 65,
    item := "Theorem 3.4, (iii) (F2) — pre-step の場合",
    sectionId := "frdi-thm-3-4" }

/-! ## ★★(F2) の Frobenius 型の成分

★★`prop_1_10_i_exists_frobType`(在庫)は pre-step 版と**同じ形の四角形**を与えるが、
`degFr φ' = degFr φ` を主張に含まない。
★しかし `degFr β = degFr α` と四角形から、**`ℕ+` の消約性**で導ける。 -/

/-- ★★★**四角形の残りの辺の次数も一致する**。

`β ≫ α = α' ≫ β'` で `degFr α = degFr α'` なら `degFr β = degFr β'`。 -/
theorem degFr_of_square_other {Dd : Type u} [Category.{v} Dd] {Cc : Type u2}
    [Category.{v2} Cc] {Φ₀ : MonoidOn.{v, u, w} Dd} (P : PreFrobenioid Cc Φ₀)
    {X Y W Z : Cc} (β : X ⟶ Y) (α : Y ⟶ Z) (α' : X ⟶ W) (β' : W ⟶ Z)
    (hsq : β ≫ α = α' ≫ β') (hdeg : P.degFr α = P.degFr α') :
    P.degFr β = P.degFr β' := by
  have h := congrArg P.degFr hsq
  rw [P.degFr_comp, P.degFr_comp, hdeg, mul_comm (P.degFr β') (P.degFr α')] at h
  exact mul_left_cancel h

def degFr_of_square_other.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 65,
    item := "Theorem 3.4, (iii) (F2) — 四角形の残りの辺の次数",
    sectionId := "frdi-thm-3-4" }

/-- ★★★★**(F2) の Frobenius 型の場合**。

★pre-step 版と同じ手だが、`hdeg` を `degFr_of_square_other` で作る。 -/
theorem admissible_of_frobType (Ψ : C₁ ⥤ C₂) (F₁ : FrobenioidCore P₁)
    (hFT : ∀ {X Y : C₁} (f : X ⟶ Y),
      IsFrobeniusType P₁ f → IsFrobeniusType P₂ (Ψ.map f))
    (hdegmap : ∀ {X Y Z W : C₁} (f : X ⟶ Y) (g : Z ⟶ W),
      P₁.degFr f = P₁.degFr g → P₂.degFr (Ψ.map f) = P₂.degFr (Ψ.map g))
    {p₁ p₂ : ℕ+} {A B : C₁} (ζ : A ⟶ B) (hζ : IsFrobeniusType P₁ ζ)
    (hA : IsAdmissibleObj P₁ P₂ Ψ p₁ p₂ A)
    {A' : C₁} (α : A ⟶ A') (hα : IsFrobeniusType P₁ α) (hdα : P₁.degFr α = p₁) :
    IsAdmissibleObj P₁ P₂ Ψ p₁ p₂ B := by
  obtain ⟨B', β, φ', hβ, hdβ, hφ', hsq⟩ :=
    prop_1_10_i_exists_frobType P₁ F₁ ζ hζ α hα
  intro E θ hθ hdθ
  refine ⟨hFT θ hθ, ?_⟩
  have h1 : P₂.degFr (Ψ.map β) = P₂.degFr (Ψ.map θ) :=
    degFr_map_eq_of_sameDeg F₁ Ψ β θ hβ hθ (by rw [hdβ, hdα, hdθ])
  have hdφ' : P₁.degFr ζ = P₁.degFr φ' :=
    degFr_of_square_other P₁ ζ β α φ' hsq hdβ
  have h2 : P₂.degFr (Ψ.map β) = P₂.degFr (Ψ.map α) := by
    have hsq' : Ψ.map ζ ≫ Ψ.map β = Ψ.map α ≫ Ψ.map φ' := by
      rw [← Ψ.map_comp, ← Ψ.map_comp, hsq]
    exact prop_1_14_iv_degFr P₂ (Ψ.map ζ) (Ψ.map β) (Ψ.map α) (Ψ.map φ') hsq'
      (hdegmap ζ φ' hdφ')
  rw [← h1, h2]
  exact (hA α hα hdα).2

def admissible_of_frobType.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 65,
    item := "Theorem 3.4, (iii) (F2) — Frobenius 型の場合",
    sectionId := "frdi-thm-3-4" }

/-- ★★★★**(F2) の伝搬の逆向き** —— 四角形は対称だから同じ手で出る。

★原典の (F2) は「**⇔**」なので、両向きが必要である。
★`arbFactor` で束ねるとき pull-back の成分だけ向きが逆になるので、これで埋める。 -/
theorem admissible_of_square' (Ψ : C₁ ⥤ C₂) (F₁ : FrobenioidCore P₁)
    {p₁ p₂ : ℕ+} {A B : C₁} (ζ : B ⟶ A) (φ : A ⟶ A) (ψ : B ⟶ B)
    (hsq : ζ ≫ φ = ψ ≫ ζ)
    (hφ : IsFrobeniusType P₁ φ) (hdφ : P₁.degFr φ = p₁)
    (hψ : IsFrobeniusType P₁ ψ) (hdψ : P₁.degFr ψ = p₁)
    (hFT : ∀ {X Y : C₁} (f : X ⟶ Y),
      IsFrobeniusType P₁ f → IsFrobeniusType P₂ (Ψ.map f))
    (hB : IsAdmissibleObj P₁ P₂ Ψ p₁ p₂ B) :
    IsAdmissibleObj P₁ P₂ Ψ p₁ p₂ A := by
  intro E θ hθ hdθ
  refine ⟨hFT θ hθ, ?_⟩
  have h1 : P₂.degFr (Ψ.map φ) = P₂.degFr (Ψ.map θ) :=
    degFr_map_eq_of_sameDeg F₁ Ψ φ θ hφ hθ (by rw [hdφ, hdθ])
  have h2 : P₂.degFr (Ψ.map φ) = P₂.degFr (Ψ.map ψ) :=
    degFr_map_eq_of_square Ψ ζ φ ψ hsq
  have h3 : P₂.degFr (Ψ.map ψ) = p₂ := (hB ψ hψ hdψ).2
  rw [← h1, h2, h3]

def admissible_of_square'.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 65,
    item := "Theorem 3.4, (iii) (F2) — 伝搬の逆向き",
    sectionId := "frdi-thm-3-4" }

/-! ## ★★★★★(F2) 全体 —— 場合分けは不要だった

★★★**測定(2026-08-18)**: 在庫の `prop_1_10_i_exists`(`Prop110.lean:424`)は
**任意の射**に対して四角形を与える(中で `arbFactor` の 3 成分を順に使っている)。
★したがって (F2) は**場合分けずに 1 本で書ける**。
★先に取った 3 成分の補題は、測定の記録として残す。 -/

/-- ★★★★★**[FrdI] Theorem 3.4, (iii) の (F2)** ——
**admissible 性は任意の射を越えて伝わる**。 -/
theorem admissible_of_hom (Ψ : C₁ ⥤ C₂) (F₁ : FrobenioidCore P₁)
    (hFT : ∀ {X Y : C₁} (f : X ⟶ Y),
      IsFrobeniusType P₁ f → IsFrobeniusType P₂ (Ψ.map f))
    (hdegmap : ∀ {X Y Z W : C₁} (f : X ⟶ Y) (g : Z ⟶ W),
      P₁.degFr f = P₁.degFr g → P₂.degFr (Ψ.map f) = P₂.degFr (Ψ.map g))
    {p₁ p₂ : ℕ+} {A B : C₁} (ζ : A ⟶ B)
    (hA : IsAdmissibleObj P₁ P₂ Ψ p₁ p₂ A)
    {A' : C₁} (α : A ⟶ A') (hα : IsFrobeniusType P₁ α) (hdα : P₁.degFr α = p₁) :
    IsAdmissibleObj P₁ P₂ Ψ p₁ p₂ B := by
  obtain ⟨B', β, φ', hβ, hdβ, hsq⟩ := prop_1_10_i_exists P₁ F₁ ζ α hα
  intro E θ hθ hdθ
  refine ⟨hFT θ hθ, ?_⟩
  have h1 : P₂.degFr (Ψ.map β) = P₂.degFr (Ψ.map θ) :=
    degFr_map_eq_of_sameDeg F₁ Ψ β θ hβ hθ (by rw [hdβ, hdα, hdθ])
  have hdζ : P₁.degFr ζ = P₁.degFr φ' :=
    degFr_of_square_other P₁ ζ β α φ' hsq hdβ
  have h2 : P₂.degFr (Ψ.map β) = P₂.degFr (Ψ.map α) := by
    have hsq' : Ψ.map ζ ≫ Ψ.map β = Ψ.map α ≫ Ψ.map φ' := by
      rw [← Ψ.map_comp, ← Ψ.map_comp, hsq]
    exact prop_1_14_iv_degFr P₂ (Ψ.map ζ) (Ψ.map β) (Ψ.map α) (Ψ.map φ') hsq'
      (hdegmap ζ φ' hdζ)
  rw [← h1, h2]
  exact (hA α hα hdα).2

def admissible_of_hom.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 64,
    item := "Theorem 3.4, (iii) — (F2)",
    sectionId := "frdi-thm-3-4" }

/-! ## ★★★(F3) の数論的な段 —— 順序を保つ素数の全単射は恒等

原文 (FrdI p.65):
> are prime-Frobenius base-identity endomorphisms such that degFr

★★原典は「Ψ が素数次数の**大小**を保つ」ことから `Ψ_N = id` を出す。
★その最後の段は**純粋に数論的**である。 -/

/-- ★★★★**順序を保つ素数の全単射は恒等**。

★手: 最小の反例を取る。`σ p < p` なら帰納仮定と単射性で矛盾、
`p < σ p` なら全射性で `p = σ q` を取って三分法で矛盾。 -/
theorem primeEquiv_eq_self (σ : Nat.Primes ≃ Nat.Primes)
    (hmono : ∀ {p q : Nat.Primes}, (p : ℕ) < (q : ℕ) → ((σ p : ℕ)) < ((σ q : ℕ))) :
    ∀ p : Nat.Primes, σ p = p := by
  have key : ∀ n : ℕ, ∀ p : Nat.Primes, (p : ℕ) = n → σ p = p := by
    intro n
    induction n using Nat.strong_induction_on with
    | _ n ih =>
      intro p hpn
      rcases lt_trichotomy ((σ p : ℕ)) ((p : ℕ)) with h | h | h
      · have h1 : σ (σ p) = σ p := ih ((σ p : ℕ)) (hpn ▸ h) (σ p) rfl
        have h2 : σ p = p := σ.injective h1
        exact absurd (congrArg (fun x : Nat.Primes => (x : ℕ)) h2) (ne_of_lt h)
      · exact Subtype.ext h
      · obtain ⟨q, hq⟩ := σ.surjective p
        rcases lt_trichotomy ((q : ℕ)) ((p : ℕ)) with hq1 | hq1 | hq1
        · have : σ q = q := ih ((q : ℕ)) (hpn ▸ hq1) q rfl
          rw [hq] at this
          exact absurd (congrArg (fun x : Nat.Primes => (x : ℕ)) this.symm) (ne_of_lt hq1)
        · have : q = p := Subtype.ext hq1
          rw [this] at hq
          exact absurd (congrArg (fun x : Nat.Primes => (x : ℕ)) hq) (ne_of_gt h)
        · have := hmono hq1
          rw [hq] at this
          exact absurd this (not_lt.mpr (le_of_lt h))
  exact fun p => key (p : ℕ) p rfl

def primeEquiv_eq_self.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 65,
    item := "Theorem 3.4, (iii) (F3) — 順序を保つ素数の全単射は恒等",
    sectionId := "frdi-thm-3-4" }

/-! ## ★★★(F3) の順序の中核 —— 単系の補題に落ちる

★★★**測定(2026-08-18)**: `prop_1_14_v_mp` の中で
**`Div α' = (degFr ψ) • Div α`**(`Prop114.lean:1499`)が作られている。
★原典の「step `β ∘ α`」はこの `α'` であり、その `Div` は
**次数の倍**になる。

★したがって「`β_θ = γ ∘ β_ζ` となる step `γ` がある」は
`MLe (p • a) (q • a)` に帰着し、さらに **`p ≤ q`** に帰着する。 -/

/-- ★★★★**`p • a ≤ q • a ⇔ p ≤ q`**(`a ≠ 0`、sharp かつ integral)。

★これが (F3) の「次数の大小」の中核である。 -/
theorem nsmul_mle_nsmul_iff {M : Type w} [AddCommMonoid M]
    (hint : IsIntegralMonoid M) (hsharp : IsSharp M) {a : M} (ha : a ≠ 0) (p q : ℕ) :
    MLe (p • a) (q • a) ↔ p ≤ q := by
  letI := isCancelAdd_of_isIntegralMonoid (M := M) hint
  constructor
  · rintro ⟨c, hc⟩
    by_contra hlt
    push Not at hlt
    obtain ⟨r, hr⟩ : ∃ r, p = q + r + 1 := ⟨p - q - 1, by omega⟩
    rw [hr, add_smul, add_smul, one_smul, add_assoc, add_assoc] at hc
    have hz : r • a + (a + c) = 0 := by
      have h0 : q • a + (r • a + (a + c)) = q • a + 0 := by rw [add_zero]; exact hc
      exact add_left_cancel h0
    have h1 := eq_zero_of_add_eq_zero_of_isSharp hsharp hz
    have h2 := eq_zero_of_add_eq_zero_of_isSharp hsharp h1.2
    exact ha h2.1
  · intro h
    obtain ⟨r, hr⟩ : ∃ r, q = p + r := ⟨q - p, by omega⟩
    exact ⟨r • a, by rw [hr, add_smul]⟩

def nsmul_mle_nsmul_iff.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 65,
    item := "Theorem 3.4, (iii) (F3) — 次数の大小の単系版",
    sectionId := "frdi-thm-3-4" }

/-- ★★★★**(F3) の順序の特徴づけ**。

★`prop_1_14_v_mp` が作る step の `Div` は `d • Div α`(d は Frobenius 次数)なので、
**step の分解可能性 ⇔ 次数の大小**となる。 -/
theorem step_mle_iff_degFr_le (hiso : ∀ X : C₁, IsIsotropic P₁ X)
    {A B B₁ B₂ : C₁} {p q : ℕ}
    (α : A ⟶ B) (hαs : IsStep P₁ α)
    (α₁ : A ⟶ B₁) (α₂ : A ⟶ B₂)
    (hd₁ : P₁.Div α₁ = p • P₁.Div α) (hd₂ : P₁.Div α₂ = q • P₁.Div α) :
    MLe (P₁.Div α₁) (P₁.Div α₂) ↔ p ≤ q := by
  have hdα : P₁.Div α ≠ 0 := fun h => hαs.2 (hiso A B α h hαs.1)
  rw [hd₁, hd₂]
  exact nsmul_mle_nsmul_iff (P₁.divisorial (P₁.toElem.obj A).base).1.1
    (P₁.divisorial (P₁.toElem.obj A).base).2 hdα p q

def step_mle_iff_degFr_le.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 65,
    item := "Theorem 3.4, (iii) (F3) — step の分解可能性 ⇔ 次数の大小",
    sectionId := "frdi-thm-3-4" }

/-! ## ★★★(F3) の順序の保存

★★`step_mle_iff_degFr_le` で「次数の大小」が
**step の分解可能性**として圧倒的に書けたので、
`Ψ` が step を保てば順序も保たれる。 -/

/-- ★★★★**step の分解可能性は関手で移る**。

★これが (F3) の「Ψ が順序を保つ」の中核である。 -/
theorem stepDecomp_map (Ψ : C₁ ⥤ C₂)
    (hstep : ∀ {X Y : C₁} (f : X ⟶ Y), IsStep P₁ f → IsStep P₂ (Ψ.map f))
    {A B₁ B₂ : C₁} (α₁ : A ⟶ B₁) (α₂ : A ⟶ B₂)
    (h : ∃ γ : B₁ ⟶ B₂, IsStep P₁ γ ∧ α₁ ≫ γ = α₂) :
    ∃ γ : Ψ.obj B₁ ⟶ Ψ.obj B₂, IsStep P₂ γ ∧ Ψ.map α₁ ≫ γ = Ψ.map α₂ := by
  obtain ⟨γ, hγ, hcomp⟩ := h
  exact ⟨Ψ.map γ, hstep γ hγ, by rw [← Ψ.map_comp, hcomp]⟩

def stepDecomp_map.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 65,
    item := "Theorem 3.4, (iii) (F3) — step の分解可能性の移送",
    sectionId := "frdi-thm-3-4" }

/-- ★★★★★**(F3) の順序の保存**。

★両側の `step_mle_iff_degFr_le` を繋ぐ。
★橋渡し(`hbridge`)は `stepDecomp_map` が与える。 -/
theorem degFr_order_preserve (Ψ : C₁ ⥤ C₂)
    (hiso₁ : ∀ X : C₁, IsIsotropic P₁ X) (hiso₂ : ∀ X : C₂, IsIsotropic P₂ X)
    {A B B₁ B₂ : C₁} {p q p' q' : ℕ}
    (α : A ⟶ B) (hαs : IsStep P₁ α) (hΨα : IsStep P₂ (Ψ.map α))
    (α₁ : A ⟶ B₁) (α₂ : A ⟶ B₂)
    (hd₁ : P₁.Div α₁ = p • P₁.Div α) (hd₂ : P₁.Div α₂ = q • P₁.Div α)
    (hd₁' : P₂.Div (Ψ.map α₁) = p' • P₂.Div (Ψ.map α))
    (hd₂' : P₂.Div (Ψ.map α₂) = q' • P₂.Div (Ψ.map α))
    (hbridge : MLe (P₁.Div α₁) (P₁.Div α₂) ↔
      MLe (P₂.Div (Ψ.map α₁)) (P₂.Div (Ψ.map α₂))) :
    p ≤ q ↔ p' ≤ q' := by
  rw [← step_mle_iff_degFr_le hiso₁ α hαs α₁ α₂ hd₁ hd₂,
    ← step_mle_iff_degFr_le hiso₂ (Ψ.map α) hΨα (Ψ.map α₁) (Ψ.map α₂) hd₁' hd₂']
  exact hbridge

def degFr_order_preserve.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 65,
    item := "Theorem 3.4, (iii) (F3) — 順序の保存",
    sectionId := "frdi-thm-3-4" }

/-- ★★★★**(F2) の逆向き(任意の射)** —— zigzag を渡るのに要る。

★同じ四角形を使い、`hA` を `hB` に入れ替えるだけ。 -/
theorem admissible_of_hom' (Ψ : C₁ ⥤ C₂) (F₁ : FrobenioidCore P₁)
    (hFT : ∀ {X Y : C₁} (f : X ⟶ Y),
      IsFrobeniusType P₁ f → IsFrobeniusType P₂ (Ψ.map f))
    (hdegmap : ∀ {X Y Z W : C₁} (f : X ⟶ Y) (g : Z ⟶ W),
      P₁.degFr f = P₁.degFr g → P₂.degFr (Ψ.map f) = P₂.degFr (Ψ.map g))
    {p₁ p₂ : ℕ+} {A B : C₁} (ζ : A ⟶ B)
    (hB : IsAdmissibleObj P₁ P₂ Ψ p₁ p₂ B)
    {A' : C₁} (α : A ⟶ A') (hα : IsFrobeniusType P₁ α) (hdα : P₁.degFr α = p₁) :
    IsAdmissibleObj P₁ P₂ Ψ p₁ p₂ A := by
  obtain ⟨B', β, φ', hβ, hdβ, hsq⟩ := prop_1_10_i_exists P₁ F₁ ζ α hα
  intro E θ hθ hdθ
  refine ⟨hFT θ hθ, ?_⟩
  have h1 : P₂.degFr (Ψ.map α) = P₂.degFr (Ψ.map θ) :=
    degFr_map_eq_of_sameDeg F₁ Ψ α θ hα hθ (by rw [hdα, hdθ])
  have hdζ : P₁.degFr ζ = P₁.degFr φ' :=
    degFr_of_square_other P₁ ζ β α φ' hsq hdβ
  have h2 : P₂.degFr (Ψ.map β) = P₂.degFr (Ψ.map α) := by
    have hsq' : Ψ.map ζ ≫ Ψ.map β = Ψ.map α ≫ Ψ.map φ' := by
      rw [← Ψ.map_comp, ← Ψ.map_comp, hsq]
    exact prop_1_14_iv_degFr P₂ (Ψ.map ζ) (Ψ.map β) (Ψ.map α) (Ψ.map φ') hsq'
      (hdegmap ζ φ' hdζ)
  have h3 : P₂.degFr (Ψ.map β) = p₂ :=
    (hB β hβ (by rw [hdβ, hdα])).2
  rw [← h1, ← h2, h3]

def admissible_of_hom'.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 64,
    item := "Theorem 3.4, (iii) (F2) — 逆向き",
    sectionId := "frdi-thm-3-4" }

/-- ★★★★★**[FrdI] Theorem 3.4, (iii) の (F2)** —— 原典どおりの**同値**の形。

原文 (FrdI p.64):
> B1 in C1, A1 is (p1, p2)-admissible if and only if B1 is. -/
theorem admissible_iff_of_hom (Ψ : C₁ ⥤ C₂) (F₁ : FrobenioidCore P₁)
    (hFT : ∀ {X Y : C₁} (f : X ⟶ Y),
      IsFrobeniusType P₁ f → IsFrobeniusType P₂ (Ψ.map f))
    (hdegmap : ∀ {X Y Z W : C₁} (f : X ⟶ Y) (g : Z ⟶ W),
      P₁.degFr f = P₁.degFr g → P₂.degFr (Ψ.map f) = P₂.degFr (Ψ.map g))
    {p₁ p₂ : ℕ+} {A B : C₁} (ζ : A ⟶ B)
    {A' : C₁} (α : A ⟶ A') (hα : IsFrobeniusType P₁ α) (hdα : P₁.degFr α = p₁) :
    IsAdmissibleObj P₁ P₂ Ψ p₁ p₂ A ↔ IsAdmissibleObj P₁ P₂ Ψ p₁ p₂ B :=
  ⟨fun hA => admissible_of_hom Ψ F₁ hFT hdegmap ζ hA α hα hdα,
   fun hB => admissible_of_hom' Ψ F₁ hFT hdegmap ζ hB α hα hdα⟩

def admissible_iff_of_hom.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 64,
    item := "Theorem 3.4, (iii) — (F2) の同値の形",
    sectionId := "frdi-thm-3-4" }

/-! ## ★★★★連結性で全対象へ広げる

原文 (FrdI p.65):
> Observe, moreover, that since C1 is connected, to complete the proof of assertion

★★両向きの (F2)(`admissible_iff_of_hom`)が揃ったので、
mathlib の `induct_on_objects` に乗せるだけ。 -/

/-- ★★★★★**連結性で admissible 性が全対象へ広がる**。 -/
theorem admissible_all_of_connected (Ψ : C₁ ⥤ C₂) (F₁ : FrobenioidCore P₁)
    (hFT : ∀ {X Y : C₁} (f : X ⟶ Y),
      IsFrobeniusType P₁ f → IsFrobeniusType P₂ (Ψ.map f))
    (hdegmap : ∀ {X Y Z W : C₁} (f : X ⟶ Y) (g : Z ⟶ W),
      P₁.degFr f = P₁.degFr g → P₂.degFr (Ψ.map f) = P₂.degFr (Ψ.map g))
    {p₁ p₂ : ℕ+} {A₀ : C₁} (hA₀ : IsAdmissibleObj P₁ P₂ Ψ p₁ p₂ A₀)
    (hex : ∀ X : C₁, ∃ (X' : C₁) (α : X ⟶ X'),
      IsFrobeniusType P₁ α ∧ P₁.degFr α = p₁)
    (A : C₁) : IsAdmissibleObj P₁ P₂ Ψ p₁ p₂ A := by
  haveI := P₁.connectedC
  have key : ∀ X : C₁, X ∈ {X : C₁ | IsAdmissibleObj P₁ P₂ Ψ p₁ p₂ X} := by
    intro X
    refine CategoryTheory.induct_on_objects {X : C₁ | IsAdmissibleObj P₁ P₂ Ψ p₁ p₂ X}
      (show A₀ ∈ {X : C₁ | IsAdmissibleObj P₁ P₂ Ψ p₁ p₂ X} from hA₀) ?_ X
    intro j₁ j₂ f
    obtain ⟨X', α, hα, hdα⟩ := hex j₁
    exact admissible_iff_of_hom Ψ F₁ hFT hdegmap f α hα hdα
  intro B φ hφ hdφ
  exact key A φ hφ hdφ

def admissible_all_of_connected.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 65,
    item := "Theorem 3.4, (iii) — 連結性で全対象へ",
    sectionId := "frdi-thm-3-4" }

/-! ## ★★素数の対応の well-defined 性

★(F1) は「ある `p₂`」を与えるだけなので、
**それが一意である**ことを別に言う必要がある。 -/

/-- ★★★★**admissible の `p₂` は一意**。

★同じ対象が `(p₁,p₂)`と`(p₁,p₂')` の両方で admissible なら `p₂ = p₂'`。 -/
theorem admissible_snd_unique (Ψ : C₁ ⥤ C₂) {p₁ p₂ p₂' : ℕ+} {A : C₁}
    (h : IsAdmissibleObj P₁ P₂ Ψ p₁ p₂ A) (h' : IsAdmissibleObj P₁ P₂ Ψ p₁ p₂' A)
    {A' : C₁} (α : A ⟶ A') (hα : IsFrobeniusType P₁ α) (hdα : P₁.degFr α = p₁) :
    p₂ = p₂' := by
  rw [← (h α hα hdα).2, ← (h' α hα hdα).2]

def admissible_snd_unique.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 64,
    item := "Theorem 3.4, (iii) — 素数の対応の一意性",
    sectionId := "frdi-thm-3-4" }

/-! ## ★★次数の対応を任意の次数へ広げる

原文 (FrdI p.64):
> degree ΨN≥1 (d) [i.e., since arbitrary morphisms may be written as composites of

★★原典は「任意の射は prime-Frobenius と linear の合成で書ける」ことを使う。
★在庫の `prop_1_10_v`(Frobenius 型 ⇔ prime-Frobenius の合成)がそれを与える。
★合成で次数が掛け算になることは関手性と `degFr_comp` から直ちに出る。 -/

/-- ★★★**像の次数は合成で掛け算になる**。 -/
theorem degFr_map_comp (Ψ : C₁ ⥤ C₂) {A B E : C₁} (φ : A ⟶ B) (ψ : B ⟶ E) :
    P₂.degFr (Ψ.map (φ ≫ ψ)) = P₂.degFr (Ψ.map ψ) * P₂.degFr (Ψ.map φ) := by
  rw [Ψ.map_comp, P₂.degFr_comp]

def degFr_map_comp.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 64,
    item := "Theorem 3.4, (iii) — 像の次数の乗法性",
    sectionId := "frdi-thm-3-4" }

/-- ★★★★★**[FrdI] Theorem 3.4, (iii)** の次数の対応 ——
**Ψ は次数 d の射を次数 Ψ_N(d) の射へ送る**。

原文 (FrdI p.64):
> degree ΨN≥1 (d) [i.e., since arbitrary morphisms may be written as composites of

★`IsPrimeFrobComposite` の帰納法。`iso` の枝は次数 1、
`cons` の枝は `degFr_map_comp` と `Ψ_N` の乗法性。 -/
theorem degFr_map_of_frobComposite (Ψ : C₁ ⥤ C₂) (ΨN : ℕ+ →* ℕ+)
    (hisomap : ∀ {X Y : C₁} (f : X ⟶ Y), IsIso f → IsIso (Ψ.map f))
    (hprime : ∀ {X Y : C₁} (f : X ⟶ Y), IsPrimeFrobenius P₁ f →
      P₂.degFr (Ψ.map f) = ΨN (P₁.degFr f))
    {A B : C₁} {φ : A ⟶ B} (h : IsPrimeFrobComposite P₁ φ) :
    P₂.degFr (Ψ.map φ) = ΨN (P₁.degFr φ) := by
  induction h with
  | iso f hf =>
    haveI := hf
    haveI := hisomap f hf
    rw [degFr_of_isIso P₂ (Ψ.map f), degFr_of_isIso P₁ f, map_one]
  | cons hφ _ ih =>
    rw [degFr_map_comp, P₁.degFr_comp, map_mul, ih, hprime _ hφ]

def degFr_map_of_frobComposite.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 64,
    item := "Theorem 3.4, (iii) — 次数 d を Ψ_N(d) へ送る",
    sectionId := "frdi-thm-3-4" }

/-! ## ★★★★`Ψ_N = id` —— (F3) の結論

原文 (FrdI p.62):
> if C1, C2 admit a non-group-like object, then ΨN≥1 is the identity -/

/-- ★★★★★**素数で恒等なら `Ψ_N` も恒等**。

★`F3-identity`(順序を保つ素数の全単射は恒等)と
`psiN-extend`(素数の全単射 ⇒ `ℕ≥1 ≃* ℕ≥1`)を繋ぐ。 -/
theorem pnatMulEquivOfPrimeEquiv_eq_self (σ : Nat.Primes ≃ Nat.Primes)
    (hσ : ∀ p, σ p = p) (n : ℕ+) : pnatMulEquivOfPrimeEquiv σ n = n := by
  show PrimeMultiset.prod (Multiset.map σ (PNat.factorMultiset n)) = n
  rw [show (σ : Nat.Primes → Nat.Primes) = id from funext hσ, Multiset.map_id]
  exact PNat.prod_factorMultiset n

/-- ★★★★★**[FrdI] Theorem 3.4, (iii)** の「非 group-like なら `Ψ_N` は恒等」。

★順序を保つことだけから出る。 -/
theorem psiN_eq_id_of_orderPreserve (σ : Nat.Primes ≃ Nat.Primes)
    (hmono : ∀ {p q : Nat.Primes}, (p : ℕ) < (q : ℕ) → ((σ p : ℕ)) < ((σ q : ℕ)))
    (n : ℕ+) : pnatMulEquivOfPrimeEquiv σ n = n :=
  pnatMulEquivOfPrimeEquiv_eq_self σ (primeEquiv_eq_self σ hmono) n

def psiN_eq_id_of_orderPreserve.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 62,
    item := "Theorem 3.4, (iii) — 非 group-like なら Ψ_N は恒等",
    sectionId := "frdi-thm-3-4" }

/-! ## ★★★素数の対応の全単射性 —— quasi-inverse で往復する

原文 (FrdI p.64):
> which [by considering a quasi-inverse to Ψ] is easily seen to be an automorphism.

★★`Ψ` で `p₁ ↦ p₂`、`Ψ⁻¹` で `p₂ ↦ p₃` とすると `p₃ = p₁` になる。 -/

/-- ★★★★**往復すると元に戻る**。

★これが素数の対応が**全単射**であることの中核である。 -/
theorem admissible_roundtrip (e : C₁ ≌ C₂) {p₁ p₂ p₃ : ℕ+} {A : C₁}
    (h₁ : IsAdmissibleObj P₁ P₂ e.functor p₁ p₂ A)
    (h₂ : IsAdmissibleObj P₂ P₁ e.inverse p₂ p₃ (e.functor.obj A))
    {A' : C₁} (α : A ⟶ A') (hα : IsFrobeniusType P₁ α) (hdα : P₁.degFr α = p₁)
    (hdeg : P₁.degFr (e.inverse.map (e.functor.map α)) = P₁.degFr α) :
    p₃ = p₁ := by
  have hb := h₁ α hα hdα
  have h := (h₂ (e.functor.map α) hb.1 hb.2).2
  rw [hdeg, hdα] at h
  exact h.symm

def admissible_roundtrip.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 64,
    item := "Theorem 3.4, (iii) — quasi-inverse で往復すると元に戻る",
    sectionId := "frdi-thm-3-4" }

/-! ## ★★★★射クラスの保存 —— まず Frobenius 型

原文 (FrdI p.64):
> tion 1.10, (v), the condition of the claim implies that Ψistr

★★prime-Frobenius の保存から、`prop_1_10_v`
(Frobenius 型 ⇔ prime-Frobenius の合成)で Frobenius 型の保存が出る。 -/

/-- ★★★**prime-Frobenius の合成は関手で保たれる**。 -/
theorem isPrimeFrobComposite_map (Ψ : C₁ ⥤ C₂)
    (hisomap : ∀ {X Y : C₁} (f : X ⟶ Y), IsIso f → IsIso (Ψ.map f))
    (hprime : ∀ {X Y : C₁} (f : X ⟶ Y),
      IsPrimeFrobenius P₁ f → IsPrimeFrobenius P₂ (Ψ.map f))
    {A B : C₁} {φ : A ⟶ B} (h : IsPrimeFrobComposite P₁ φ) :
    IsPrimeFrobComposite P₂ (Ψ.map φ) := by
  induction h with
  | iso f hf => exact IsPrimeFrobComposite.iso _ (hisomap f hf)
  | cons hφ _ ih =>
    rw [Ψ.map_comp]
    exact IsPrimeFrobComposite.cons (hprime _ hφ) ih

/-- ★★★★**Ψ は Frobenius 型を保つ**。

★`prop_1_10_v` を両側で使う。 -/
theorem isFrobeniusType_map_of_primeFrob (Ψ : C₁ ⥤ C₂)
    (F₁ : FrobenioidCore P₁) (F₂ : FrobenioidCore P₂)
    (hisomap : ∀ {X Y : C₁} (f : X ⟶ Y), IsIso f → IsIso (Ψ.map f))
    (hprime : ∀ {X Y : C₁} (f : X ⟶ Y),
      IsPrimeFrobenius P₁ f → IsPrimeFrobenius P₂ (Ψ.map f))
    {A B : C₁} (φ : A ⟶ B) (h : IsFrobeniusType P₁ φ) :
    IsFrobeniusType P₂ (Ψ.map φ) :=
  (prop_1_10_v P₂ F₂ (Ψ.map φ)).mpr
    (isPrimeFrobComposite_map Ψ hisomap hprime ((prop_1_10_v P₁ F₁ φ).mp h))

def isFrobeniusType_map_of_primeFrob.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 64,
    item := "Theorem 3.4, (iii) — Frobenius 型の保存",
    sectionId := "frdi-thm-3-4" }

end ABC3.Found.FrdI
