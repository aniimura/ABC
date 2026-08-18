/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Prop32Equiv

/-!
# [FrdI] `Proposition 3.2, (ii)` —— `𝒞^pf` の射の辞書(10 種)

原文 (FrdI p.59):
> (ii) An arrow of Cpf is a(n) morphism of Frobenius type (respectively,

原文 (FrdI p.59):
> pre-step; base-isomorphism; base-identity endomorphism; isomorphism;

原文 (FrdI p.59):
> pull-back morphism; isometry; co-angular morphism; LB-invertible mor-

原文 (FrdI p.59):
> phism; morphism of a given Frobenius degree) if and only if a cofinal collection

## ★★「共終な collection」の読み方(実装上の判断)

★原文は「**共終な**代表元の族がそうであること」と言う。`𝒞^pf` の射は
添字圏(filtered)上の余極限なので、これは
「**添字を十分押し上げた代表元**でそうであること」である。

★★`𝒞` が Frobenius-isotropic 型なら、添字は**始域が isotropic な所まで押し上げられる**
(`exists_idx_isotropic` → `exists_rep_isotropic`)。★そこでは

- **co-angular が自動**(`isCoAngular_of_isotropic_dom`)、
- **等長な pre-step は同型**(`IsIsotropic` の定義)

なので、10 種の判定はすべて**代表元のそのままの判定と一致する**。
★したがって本ファイルは各条を**2 段**で置く:

1. `…_mk_iff` —— 任意の代表元での判定(`𝒞^pf` 側で自動になる条を消した形)
2. `…_mk_iff_of_isotropic` —— **共終な**(始域 isotropic の)代表元での判定。
   ここで初めて「`𝒞` の同じ述語」と一致する。

## ★10 種の対応表

| 原文の語 | 宣言 | 置き場所 |
|---|---|---|
| morphism of Frobenius type | `isFrobeniusType_mk_iff(_of_isotropic)` | 本ファイル |
| pre-step | `isPreStep_mk_iff` | `Prop32.lean` |
| base-isomorphism | `isBaseIsomorphism_mk_iff` | `Prop32.lean` |
| base-identity endomorphism | ★`isBaseIdentity_mk_diag` ＋ `exists_rep_diag` | 本ファイル |
| isomorphism | `pfRoot_isIso_mk_iff(_of_isotropic)` | 本ファイル |
| pull-back morphism | `isPullBack_mk_iff(_of_isotropic)` | 本ファイル |
| isometry | `isIsometric_mk_iff` | `Prop32.lean` |
| co-angular morphism | `isCoAngular_mk_iff_of_isotropic` ＋ `exists_rep_coAngular` | 本ファイル |
| LB-invertible morphism | `isLBInvertible_mk_iff(_of_isotropic)` | 本ファイル |
| morphism of a given Frobenius degree | `degFr_mk_iff` | `Prop32.lean` |

## ★★測定 —— pull-back は 3 つの述語に分解できる

★`Definition 1.3, (iv)(b)`(`pullBackLB`)と `Proposition 1.4, (ii)`
(`prop_1_4_ii_mpr`)を合わせると、**Frobenioid では**

  `IsPullBack φ ⟺ co-angular ∧ 等長 ∧ linear`

である(`isPullBack_iff`)。★これで pull-back の条も他の 9 条と同じ道具で落ちる。
★普遍性(`Function.Bijective`)を直接扱わずに済むのが要点。

## ★base-identity だけは「対角の添字」が要る

★自己射 `f : X ⟶ X` の代表元 `φ : Z.1 ⟶ Z.2` は、`Z.1 ≠ Z.2` だと
`IsBaseIdentity P φ` が**型として書けない**。★そこで添字を
「2 脚が同じ射」の所まで押し上げる(`exists_idx_diag`)。
そこでは `repBase = (Base l) ∘ Base a ∘ (Base l)⁻¹` が共役なので、
`Base a = 𝟙` と `rootBase = 𝟙` が同値になる。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w v2 u2

/-! ## ★1. pull-back の分解(Frobenioid 一般) -/

section General

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D}

/-- ★★★**Frobenioid では pull-back は「LB-invertible ＋ linear」に他ならない**。

★→ は `Definition 1.3, (iv)(b)`(`pullBackLB`)、← は `Proposition 1.4, (ii)`。
★★これで `Proposition 3.2, (ii)` の pull-back の条が、
co-angular・等長・linear の 3 つに分解して扱えるようになる。 -/
theorem isPullBack_iff (Q : PreFrobenioid C Φ) (Fc : FrobenioidCore Q) {A B : C} (φ : A ⟶ B) :
    IsPullBack Q φ ↔ IsLBInvertible Q φ ∧ IsLinear Q φ :=
  ⟨Fc.pullBackLB φ, fun h => prop_1_4_ii_mpr Q Fc φ h.1 h.2⟩

variable {P : PreFrobenioid C Φ}

/-! ## ★2. isotropic な始域では co-angular が自動になり、判定が縮む -/

theorem isLBInvertible_iff_of_isotropic (Fc : FrobenioidCore P) {A B : C}
    (hA : IsIsotropic P A) (φ : A ⟶ B) : IsLBInvertible P φ ↔ IsIsometric P φ :=
  ⟨fun h => h.2, fun h => ⟨isCoAngular_of_isotropic_dom (P := P) Fc hA φ, h⟩⟩

theorem isFrobeniusType_iff_of_isotropic (Fc : FrobenioidCore P) {A B : C}
    (hA : IsIsotropic P A) (φ : A ⟶ B) :
    IsFrobeniusType P φ ↔ IsIsometric P φ ∧ IsBaseIsomorphism P φ :=
  and_congr (isLBInvertible_iff_of_isotropic Fc hA φ) Iff.rfl

theorem isPullBack_iff_of_isotropic (Fc : FrobenioidCore P) {A B : C}
    (hA : IsIsotropic P A) (φ : A ⟶ B) :
    IsPullBack P φ ↔ IsIsometric P φ ∧ IsLinear P φ :=
  (isPullBack_iff P Fc φ).trans
    (and_congr (isLBInvertible_iff_of_isotropic Fc hA φ) Iff.rfl)

theorem isIso_iff_of_isotropic {A B : C} (hA : IsIsotropic P A) (φ : A ⟶ B) :
    IsIso φ ↔ IsPreStep P φ ∧ IsIsometric P φ :=
  ⟨fun _ => ⟨isPreStep_of_isIso P φ, isIsometric_of_isIso P φ⟩,
    fun h => hA _ φ h.2 h.1⟩

end General

section Dict

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} {P : PreFrobenioid C Φ} {F : FrobenioidCore P}

/-! ## ★3. 共終性 —— 代表元は「始域が isotropic」まで押し上げられる -/

/-- ★★**共終性** —— どの射も「始域が isotropic な代表元」を持つ。 -/
theorem exists_rep_isotropic (hfi : IsOfFrobeniusIsotropicType P) {X Y : PfRootObj P F}
    (f : X ⟶ Y) :
    ∃ (Z : IdxPf P F (rtObj P F X.obj Y.root) (rtObj P F Y.obj X.root))
      (φ : Z.right.obj.1 ⟶ Z.right.obj.2),
      f = HomPf.mk Z φ ∧ IsIsotropic P Z.right.obj.1 := by
  obtain ⟨Z₀, φ₀, rfl⟩ := HomPf.exists_rep f
  obtain ⟨Z, u, hZ⟩ := exists_idx_isotropic (F := F) hfi Z₀
  exact ⟨Z, idxTransport P F u φ₀, (HomPf.mk_map u φ₀).symm, hZ⟩

/-- ★★**(ii) co-angular の `𝒞` 側** —— 共終的に co-angular な代表元が取れる。

★`𝒞^pf` 側は `pfRoot_isCoAngular` で**常に**成り立つ。 -/
theorem exists_rep_coAngular (hfi : IsOfFrobeniusIsotropicType P) {X Y : PfRootObj P F}
    (f : X ⟶ Y) :
    ∃ (Z : IdxPf P F (rtObj P F X.obj Y.root) (rtObj P F Y.obj X.root))
      (φ : Z.right.obj.1 ⟶ Z.right.obj.2),
      f = HomPf.mk Z φ ∧ IsCoAngular P φ := by
  obtain ⟨Z, φ, hf, hiso⟩ := exists_rep_isotropic hfi f
  exact ⟨Z, φ, hf, isCoAngular_of_isotropic_dom (P := P) F hiso φ⟩

/-! ## ★4. 任意の代表元での判定(`𝒞^pf` 側で自動になる条を消した形) -/

/-- ★★**(ii) LB-invertible** —— `𝒞^pf` では co-angular が自動なので等長性だけ。 -/
theorem isLBInvertible_mk_iff (hfi : IsOfFrobeniusIsotropicType P) {X Y : PfRootObj P F}
    (Z : IdxPf P F (rtObj P F X.obj Y.root) (rtObj P F Y.obj X.root))
    (φ : Z.right.obj.1 ⟶ Z.right.obj.2) :
    IsLBInvertible (pfRootPre P F) (show HomRoot P F X Y from HomPf.mk Z φ)
      ↔ IsIsometric P φ :=
  ⟨fun h => (isIsometric_mk_iff Z φ).mp h.2,
    fun h => ⟨pfRoot_isCoAngular hfi _, (isIsometric_mk_iff Z φ).mpr h⟩⟩

/-- ★★**(ii) Frobenius 型**。 -/
theorem isFrobeniusType_mk_iff (hfi : IsOfFrobeniusIsotropicType P) {X Y : PfRootObj P F}
    (Z : IdxPf P F (rtObj P F X.obj Y.root) (rtObj P F Y.obj X.root))
    (φ : Z.right.obj.1 ⟶ Z.right.obj.2) :
    IsFrobeniusType (pfRootPre P F) (show HomRoot P F X Y from HomPf.mk Z φ)
      ↔ IsIsometric P φ ∧ IsBaseIsomorphism P φ :=
  and_congr (isLBInvertible_mk_iff hfi Z φ) (isBaseIsomorphism_mk_iff Z φ)

/-- ★★**(ii) pull-back**。 -/
theorem isPullBack_mk_iff (hfi : IsOfFrobeniusIsotropicType P) {X Y : PfRootObj P F}
    (Z : IdxPf P F (rtObj P F X.obj Y.root) (rtObj P F Y.obj X.root))
    (φ : Z.right.obj.1 ⟶ Z.right.obj.2) :
    IsPullBack (pfRootPre P F) (show HomRoot P F X Y from HomPf.mk Z φ)
      ↔ IsIsometric P φ ∧ IsLinear P φ :=
  (isPullBack_iff (pfRootPre P F) (pfRootCore hfi) _).trans
    (and_congr (isLBInvertible_mk_iff hfi Z φ) (isLinear_mk_iff Z φ))

/-- ★★**(ii) 同型** —— `𝒞^pf` は isotropic 型なので「等長な pre-step」で判定できる。 -/
theorem pfRoot_isIso_mk_iff (hfi : IsOfFrobeniusIsotropicType P) {X Y : PfRootObj P F}
    (Z : IdxPf P F (rtObj P F X.obj Y.root) (rtObj P F Y.obj X.root))
    (φ : Z.right.obj.1 ⟶ Z.right.obj.2) :
    IsIso (show X ⟶ Y from HomPf.mk Z φ) ↔ IsPreStep P φ ∧ IsIsometric P φ := by
  constructor
  · intro h
    exact ⟨(isPreStep_mk_iff Z φ).mp (isPreStep_of_isIso (pfRootPre P F)
        (show X ⟶ Y from HomPf.mk Z φ)),
      (isIsometric_mk_iff Z φ).mp (isIsometric_of_isIso (pfRootPre P F)
        (show X ⟶ Y from HomPf.mk Z φ))⟩
  · intro h
    exact pfRoot_isOfIsotropicType (F := F) hfi X Y (HomPf.mk Z φ)
      ((isIsometric_mk_iff Z φ).mpr h.2) ((isPreStep_mk_iff Z φ).mpr h.1)

/-! ## ★5. 共終な代表元での辞書 —— ここで `𝒞` の同じ述語と一致する -/

theorem isFrobeniusType_mk_iff_of_isotropic (hfi : IsOfFrobeniusIsotropicType P)
    {X Y : PfRootObj P F}
    (Z : IdxPf P F (rtObj P F X.obj Y.root) (rtObj P F Y.obj X.root))
    (φ : Z.right.obj.1 ⟶ Z.right.obj.2) (hZ : IsIsotropic P Z.right.obj.1) :
    IsFrobeniusType (pfRootPre P F) (show HomRoot P F X Y from HomPf.mk Z φ)
      ↔ IsFrobeniusType P φ :=
  (isFrobeniusType_mk_iff hfi Z φ).trans (isFrobeniusType_iff_of_isotropic F hZ φ).symm

theorem isLBInvertible_mk_iff_of_isotropic (hfi : IsOfFrobeniusIsotropicType P)
    {X Y : PfRootObj P F}
    (Z : IdxPf P F (rtObj P F X.obj Y.root) (rtObj P F Y.obj X.root))
    (φ : Z.right.obj.1 ⟶ Z.right.obj.2) (hZ : IsIsotropic P Z.right.obj.1) :
    IsLBInvertible (pfRootPre P F) (show HomRoot P F X Y from HomPf.mk Z φ)
      ↔ IsLBInvertible P φ :=
  (isLBInvertible_mk_iff hfi Z φ).trans (isLBInvertible_iff_of_isotropic F hZ φ).symm

theorem isPullBack_mk_iff_of_isotropic (hfi : IsOfFrobeniusIsotropicType P)
    {X Y : PfRootObj P F}
    (Z : IdxPf P F (rtObj P F X.obj Y.root) (rtObj P F Y.obj X.root))
    (φ : Z.right.obj.1 ⟶ Z.right.obj.2) (hZ : IsIsotropic P Z.right.obj.1) :
    IsPullBack (pfRootPre P F) (show HomRoot P F X Y from HomPf.mk Z φ)
      ↔ IsPullBack P φ :=
  (isPullBack_mk_iff hfi Z φ).trans (isPullBack_iff_of_isotropic F hZ φ).symm

theorem pfRoot_isIso_mk_iff_of_isotropic (hfi : IsOfFrobeniusIsotropicType P)
    {X Y : PfRootObj P F}
    (Z : IdxPf P F (rtObj P F X.obj Y.root) (rtObj P F Y.obj X.root))
    (φ : Z.right.obj.1 ⟶ Z.right.obj.2) (hZ : IsIsotropic P Z.right.obj.1) :
    IsIso (show X ⟶ Y from HomPf.mk Z φ) ↔ IsIso φ :=
  (pfRoot_isIso_mk_iff hfi Z φ).trans (isIso_iff_of_isotropic hZ φ).symm

theorem isCoAngular_mk_iff_of_isotropic (hfi : IsOfFrobeniusIsotropicType P)
    {X Y : PfRootObj P F}
    (Z : IdxPf P F (rtObj P F X.obj Y.root) (rtObj P F Y.obj X.root))
    (φ : Z.right.obj.1 ⟶ Z.right.obj.2) (hZ : IsIsotropic P Z.right.obj.1) :
    IsCoAngular (pfRootPre P F) (show HomRoot P F X Y from HomPf.mk Z φ)
      ↔ IsCoAngular P φ :=
  ⟨fun _ => isCoAngular_of_isotropic_dom (P := P) F hZ φ,
    fun _ => pfRoot_isCoAngular hfi _⟩

/-! ## ★6. base-identity 自己射 —— 対角の添字 -/

/-- ★★**自己射の添字を「2 脚が一致し、終域が isotropic」な所まで押し上げる**。

★`exists_idx3_diag` の 2 脚版。 -/
theorem exists_idx_diag (hfi : IsOfFrobeniusIsotropicType P) {A : C} (V : IdxPf P F A A) :
    ∃ (E : C) (l : A ⟶ E) (hl : IsFrobeniusType P l),
      IsIsotropic P E ∧ Nonempty (V ⟶ idxMk (P := P) (F := F) l l hl hl rfl) := by
  obtain ⟨hv1, hv2, h12⟩ := V.hom.property
  obtain ⟨X₃, a, c, ha, hc, had, hcd, hac⟩ :=
    frob_common_upper P F V.hom.hom.1 hv1 V.hom.hom.2 hv2
  obtain ⟨Dd, δ, hδ, hDd⟩ := hfi X₃
  have hacδ : V.hom.hom.1 ≫ (a ≫ δ) = V.hom.hom.2 ≫ (c ≫ δ) := by
    rw [← Category.assoc, ← Category.assoc, hac]
    rfl
  have hdaδ : P.degFr (a ≫ δ) = P.degFr (c ≫ δ) := by
    rw [P.degFr_comp a δ, P.degFr_comp c δ, had, hcd, h12]
    rfl
  refine ⟨Dd, V.hom.hom.1 ≫ (a ≫ δ),
    IsFrobeniusType.comp P F hv1 (IsFrobeniusType.comp P F ha hδ), hDd,
    ⟨Under.homMk (show V.right ⟶ (⟨(Dd, Dd)⟩ : BiFr P F) from
      ⟨(a ≫ δ, c ≫ δ), IsFrobeniusType.comp P F ha hδ,
        IsFrobeniusType.comp P F hc hδ, hdaδ⟩)
      (WideSubcategory.hom_ext _ (Prod.ext rfl hacδ.symm))⟩⟩

/-- ★★**(ii) base-identity の共終性** —— 自己射は対角の代表元を持つ。 -/
theorem exists_rep_diag (hfi : IsOfFrobeniusIsotropicType P) {X : PfRootObj P F} (f : End X) :
    ∃ (E : C) (l : rtObj P F X.obj X.root ⟶ E) (hl : IsFrobeniusType P l) (a : E ⟶ E),
      f = HomPf.mk (idxMk (P := P) (F := F) l l hl hl rfl) a ∧ IsIsotropic P E := by
  obtain ⟨V₀, φ₀, rfl⟩ := HomPf.exists_rep (show HomRoot P F X X from f)
  obtain ⟨E, l, hl, hE, ⟨u⟩⟩ := exists_idx_diag (F := F) hfi V₀
  exact ⟨E, l, hl, idxTransport P F u φ₀, (HomPf.mk_map u φ₀).symm, hE⟩

set_option maxHeartbeats 1000000 in
/-- ★★**(ii) base-identity 自己射** —— 対角の添字での判定。

★`oTri_mk_diag` から次数の条を落とした形。第 1・第 2 脚が**同じ射** `l` なので
`repBase = (Base l) ∘ Base a ∘ (Base l)⁻¹` が共役になり、
`Base a = 𝟙` と `rootBase = 𝟙` が同値になる。 -/
theorem isBaseIdentity_mk_diag {X : PfRootObj P F} {E : C}
    (l : rtObj P F X.obj X.root ⟶ E) (hl : IsFrobeniusType P l) (a : E ⟶ E) :
    IsBaseIdentity (pfRootPre P F)
        (show End X from HomPf.mk (idxMk (P := P) (F := F) l l hl hl rfl) a)
      ↔ IsBaseIdentity P a := by
  haveI hil : IsIso (P.Base l) := hl.2
  haveI hie : IsIso (P.Base (rtExt P F X.obj X.root)) := (rtExt_frobType P F X.obj X.root).2
  have hrep : repBase (idxMk (P := P) (F := F) l l hl hl rfl) a ≫ P.Base l
      = P.Base l ≫ P.Base a :=
    repBase_spec (idxMk (P := P) (F := F) l l hl hl rfl) a
  have hroot : rootBase (show End X from HomPf.mk (idxMk (P := P) (F := F) l l hl hl rfl) a)
        ≫ P.Base (rtExt P F X.obj X.root)
      = P.Base (rtExt P F X.obj X.root)
        ≫ pfBase (HomPf.mk (idxMk (P := P) (F := F) l l hl hl rfl) a) :=
    rootBase_spec (show End X from HomPf.mk (idxMk (P := P) (F := F) l l hl hl rfl) a)
  have hpf : pfBase (HomPf.mk (idxMk (P := P) (F := F) l l hl hl rfl) a)
      = repBase (idxMk (P := P) (F := F) l l hl hl rfl) a :=
    pfBase_mk (idxMk (P := P) (F := F) l l hl hl rfl) a
  have hidB : 𝟙 ((pfRootPre P F).toElem.obj X).base ≫ P.Base (rtExt P F X.obj X.root)
      = P.Base (rtExt P F X.obj X.root) := Category.id_comp _
  have hidL : 𝟙 ((P.toElem.obj (rtObj P F X.obj X.root)).base) ≫ P.Base l
      = P.Base l := Category.id_comp _
  constructor
  · intro hb
    have h0 : rootBase (show End X from
        HomPf.mk (idxMk (P := P) (F := F) l l hl hl rfl) a)
        = (pfRootPre P F).Base (𝟙 X) := hb
    rw [(pfRootPre P F).Base_id] at h0
    rw [h0, hpf] at hroot
    have h1 : P.Base (rtExt P F X.obj X.root) ≫ 𝟙 _
        = P.Base (rtExt P F X.obj X.root)
          ≫ repBase (idxMk (P := P) (F := F) l l hl hl rfl) a :=
      (Category.comp_id _).trans (hidB.symm.trans hroot)
    have h2 : repBase (idxMk (P := P) (F := F) l l hl hl rfl) a = 𝟙 _ :=
      ((cancel_epi (P.Base (rtExt P F X.obj X.root))).mp h1).symm
    rw [h2] at hrep
    have h3 : P.Base l ≫ 𝟙 ((P.toElem.obj E).base) = P.Base l ≫ P.Base a :=
      (Category.comp_id _).trans (hidL.symm.trans hrep)
    show P.Base a = P.Base (𝟙 E)
    rw [P.Base_id]
    exact ((cancel_epi (P.Base l)).mp h3).symm
  · intro hb
    rw [show P.Base a = 𝟙 _ from hb.trans (P.Base_id E)] at hrep
    have h1 : repBase (idxMk (P := P) (F := F) l l hl hl rfl) a ≫ P.Base l
        = 𝟙 _ ≫ P.Base l := hrep.trans ((Category.comp_id _).trans hidL.symm)
    have h2 : repBase (idxMk (P := P) (F := F) l l hl hl rfl) a = 𝟙 _ :=
      (cancel_mono (P.Base l)).mp h1
    rw [hpf, h2] at hroot
    show rootBase (show End X from HomPf.mk (idxMk (P := P) (F := F) l l hl hl rfl) a)
      = (pfRootPre P F).Base (𝟙 X)
    rw [(pfRootPre P F).Base_id]
    refine (cancel_mono (P.Base (rtExt P F X.obj X.root))).mp ?_
    exact hroot.trans ((Category.comp_id _).trans hidB.symm)

end Dict

/-! ## ★出典の紐付け(条つき) -/

def isFrobeniusType_mk_iff.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 59, item := "Proposition 3.2, (ii) — 10 種の射の辞書",
    sectionId := "frdi-prop-3-2" }

def isPullBack_iff.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 24,
    item := "Definition 1.3, (iv)(b) ＋ Proposition 1.4, (ii) — pull-back の特徴づけ",
    sectionId := "frdi-def-1-3-iv" }

end ABC3.Found.FrdI
