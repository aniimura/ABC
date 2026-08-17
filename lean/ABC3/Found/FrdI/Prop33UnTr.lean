import ABC3.Found.FrdI.UnTr
import ABC3.Found.FrdI.Prop16

/-!
# [FrdI] Proposition 3.3, (iv) の本体 —— `𝒞^un-tr` の **Frobenioid 構造**

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.60。

原文 (FrdI p.60):
> which is faithful and essentially surjective; moreover, this functor determines

★★原文は「**determines a natural structure of Frobenioid on `𝒞^un-tr`**」と
1 行で言うが、中身は **`FrobenioidCore` の 21 フィールド ＋ `Frobenioid` の 2 フィールド**
である(`UnTr.lean` が `PreFrobenioid` の 6 フィールドを既に埋めている)。

## ★★このファイルの立て付け(測定、2026-08-17)

★`UnTr P` は**型として `Istr P` そのもの**で、`istrToUnTr` は
**対象について恒等**かつ **full**(`UnTr.lean`)。★したがって:

* **存在型のフィールド**は `𝒞^istr` の証人(`istr_frobenioidCore`)を
  `istrToUnTr` で送るだけで済む見込み——★**ただしその性質が商を渡ることが要る**。
* **性質が商を渡るか**は類型ごとに違う(下の表)。

## ★★★射の類型の仕分け(`UnTr.lean` の測定 ＋ 本ファイルの追加測定)

| 類型 | 商を渡るか | 根拠 |
|---|---|---|
| `degFr = n` / base-isomorphism / isometry / pre-step | ★★**`Iff.rfl`** | `unTrToElem` は代表元に `P.toElem` を当てるだけ |
| isomorphism | ★**両向き取れる** | `unTr_isIso_iff`(`𝒞^istr` が isotropic 型) |
| ★co-angular | ★**片向きは取れる**(下記) | 因子分解を送って `unTr_isIso_iff` |
| pull-back / LB-invertible / Frobenius 型 | ★co-angular に依存 | —— |

★★**原文の言い方は「if and only if it ARISES FROM such an arrow of `𝒞^istr`」**
であり、「**どの持ち上げも**その性質を持つ」ではない。★存在型のフィールドを
埋めるのに要るのは **`𝒞^istr` の性質 ⟹ `𝒞^un-tr` の性質**の向きである。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} (P : PreFrobenioid C Φ)

/-! ## ★1. `𝒞^istr` の性質を `𝒞^un-tr` へ送る(易しい 4 類型)

★★`unTr_isPreStep_iff` などが `Iff.rfl` なので、そのまま使えばよい。
★ここでは `istrPre` 側の言い方との橋渡しだけを置く。 -/

variable {P} in
/-- ★`𝒞^istr` の射の `Base` は `𝒞` のもの。 -/
theorem istrPre_Base (Fc : FrobenioidCore P) {A B : Istr P} (f : A ⟶ B) :
    (istrPre P Fc).Base f = P.Base f.hom := rfl

variable {P} in
/-- ★`𝒞^istr` の射の `degFr` は `𝒞` のもの。 -/
theorem istrPre_degFr (Fc : FrobenioidCore P) {A B : Istr P} (f : A ⟶ B) :
    (istrPre P Fc).degFr f = P.degFr f.hom := rfl

variable {P} in
/-- ★`𝒞^istr` の射の `Div` は `𝒞` のもの。 -/
theorem istrPre_Div (Fc : FrobenioidCore P) {A B : Istr P} (f : A ⟶ B) :
    (istrPre P Fc).Div f = P.Div f.hom := rfl

/-! ## ★2. 埋まったフィールド -/

include P in
/-- ★★**(i)(b)** 底の同型は pre-step の span で持ち上がる。

★`𝒞^istr` の証人を `istrToUnTr` で送るだけ —— ★**pre-step は `Iff.rfl` で渡る**。 -/
theorem unTr_preStepSpan (Fc : FrobenioidCore P) :
    ∀ (A B : UnTr P)
      (α : ((unTrPre P Fc).toElem.obj A).base ⟶ ((unTrPre P Fc).toElem.obj B).base),
      IsIso α →
      ∃ (X : UnTr P) (φ : X ⟶ A) (ψ : X ⟶ B) (hφ : IsPreStep (unTrPre P Fc) φ),
        IsPreStep (unTrPre P Fc) ψ ∧
        α = @inv _ _ _ _ ((unTrPre P Fc).Base φ) hφ.2 ≫ (unTrPre P Fc).Base ψ := by
  intro A B α hα
  obtain ⟨X, φ, ψ, hφ, hψ, he⟩ :=
    (istr_frobenioidCore P Fc).preStepSpan (show Istr P from A) (show Istr P from B) α hα
  exact ⟨show UnTr P from X, (istrToUnTr P).map φ, (istrToUnTr P).map ψ, hφ, hψ, he⟩

include P in
/-- ★★**(vii)(b)** isotropic な対象からの射の終域も isotropic。

★★`𝒞^un-tr` の対象は `𝒞^istr` の対象、すなわち**すべて isotropic** である。 -/
theorem unTr_isotropic (Fc : FrobenioidCore P) (B : UnTr P) :
    IsIsotropic (unTrPre P Fc) B := by
  intro Dd f
  induction f using Quotient.inductionOn with
  | _ α =>
    intro hiso hstep
    haveI hα : IsIso α := (show Istr P from B).property (show Istr P from Dd).obj α hiso hstep
    haveI hmap : IsIso ((isotropicProp P).ι.map (ObjectProperty.homMk α :
        (show Istr P from B) ⟶ (show Istr P from Dd))) := hα
    haveI : IsIso (ObjectProperty.homMk α :
        (show Istr P from B) ⟶ (show Istr P from Dd)) :=
      (ObjectProperty.fullyFaithfulι (isotropicProp P)).isIso_of_isIso_map _
    let e := (istrToUnTr P).mapIso
      (asIso (ObjectProperty.homMk α : (show Istr P from B) ⟶ (show Istr P from Dd)))
    exact ⟨⟨e.inv, e.hom_inv_id, e.inv_hom_id⟩⟩

include P in
/-- ★★**(vii)(b)** isotropic な対象からの射の終域も isotropic —— ★**すべての対象が
isotropic** なので自明。 -/
theorem unTr_isotropicClosed (Fc : FrobenioidCore P) :
    ∀ {A B : UnTr P} (_ : A ⟶ B), IsIsotropic (unTrPre P Fc) A →
      IsIsotropic (unTrPre P Fc) B :=
  fun _ _ => unTr_isotropic P Fc _

include P in
/-- ★★**(vii)(a)** isotropic hull は**恒等射**でよい —— すべての対象が isotropic だから。 -/
theorem unTr_isotropicHullExists (Fc : FrobenioidCore P) :
    ∀ A : UnTr P, ∃ (B : UnTr P) (φ : A ⟶ B), IsIsotropicHull (unTrPre P Fc) φ := by
  intro A
  refine ⟨A, 𝟙 A, ?_, ?_, unTr_isotropic P Fc A, ?_⟩
  · show (unTrPre P Fc).Div (𝟙 A) = 0
    exact (unTrPre P Fc).Div_id A
  · exact isPreStep_id (unTrPre P Fc) A
  · intro Cc _ γ
    exact ⟨γ, (Category.id_comp γ).symm, fun β hβ => by rw [hβ, Category.id_comp]⟩

include P in
/-- ★★★**(vi)** 単元を除く忠実性 —— ★**`𝒞^un-tr` では自明**である。

★★`BaseEquivalent` と `MetricallyEquivalent` で `Base` と `Div` が一致し、
pre-step なので `degFr = 1` も一致する。★したがって **`𝔽_Φ` で同じ射**になり、
`𝒞^un-tr → 𝔽_Φ` が**忠実**なので 2 射は等しい。★よって `α = 1` で取れる。

★★**他所で一番苦労するフィールドが、ここでは商の定義そのものから出る。** -/
theorem unTr_faithfulUpToUnits (Fc : FrobenioidCore P) :
    ∀ {A B : UnTr P} (φ ψ : A ⟶ B), BaseEquivalent (unTrPre P Fc) φ ψ →
      MetricallyEquivalent (unTrPre P Fc) φ ψ → IsCoAngular (unTrPre P Fc) φ →
      IsPreStep (unTrPre P Fc) φ → IsCoAngular (unTrPre P Fc) ψ →
      IsPreStep (unTrPre P Fc) ψ →
      ∃ α : End B, α ∈ OTimes (unTrPre P Fc) B ∧ φ = ψ ≫ (α : B ⟶ B) := by
  intro A B φ ψ hb hm _ hsφ _ hsψ
  refine ⟨1, (OTimes (unTrPre P Fc) B).one_mem, ?_⟩
  rw [show ((1 : End B) : B ⟶ B) = 𝟙 B from rfl, Category.comp_id]
  refine (unTrToElem P).map_injective ?_
  refine ElemFrobCat.Hom.ext hb hm ?_
  show (unTrPre P Fc).degFr φ = (unTrPre P Fc).degFr ψ
  rw [show (unTrPre P Fc).degFr φ = 1 from hsφ.1, show (unTrPre P Fc).degFr ψ = 1 from hsψ.1]

include P in
/-- ★★**(v)(a)** pre-step は monomorphism。

★★`𝒞^un-tr → 𝔽_Φ` が**忠実**なので、`𝔽_Φ` の側で示して降ろす:
底は `Base` が同型ゆえ mono、次数は `ℕ+` の簡約、
零因子は `deg = 1` と `Φ(A)` の簡約性(divisorial ⟹ integral)。 -/
theorem unTr_preStepMono (Fc : FrobenioidCore P) :
    ∀ {A B : UnTr P} (φ : A ⟶ B), IsPreStep (unTrPre P Fc) φ → Mono φ := by
  intro A B φ hφ
  refine ⟨fun {Z} g h hgh => ?_⟩
  haveI hb : IsIso ((unTrPre P Fc).Base φ) := hφ.2
  letI := isCancelAdd_of_isIntegralMonoid (Φ.val ((unTrPre P Fc).toElem.obj Z).base)
    (P.divisorial _).1.1
  have e : (unTrToElem P).map g ≫ (unTrToElem P).map φ
      = (unTrToElem P).map h ≫ (unTrToElem P).map φ := by
    rw [← (unTrToElem P).map_comp, ← (unTrToElem P).map_comp, hgh]
  have hbase : (unTrPre P Fc).Base g = (unTrPre P Fc).Base h := by
    have h1 : (unTrPre P Fc).Base g ≫ (unTrPre P Fc).Base φ
        = (unTrPre P Fc).Base h ≫ (unTrPre P Fc).Base φ := congrArg ElemFrobCat.Hom.base e
    exact (cancel_mono ((unTrPre P Fc).Base φ)).mp h1
  have hdeg : (unTrPre P Fc).degFr g = (unTrPre P Fc).degFr h := by
    have h1 : (unTrPre P Fc).degFr φ * (unTrPre P Fc).degFr g
        = (unTrPre P Fc).degFr φ * (unTrPre P Fc).degFr h := congrArg ElemFrobCat.Hom.deg e
    exact mul_left_cancel h1
  have hdiv : (unTrPre P Fc).Div g = (unTrPre P Fc).Div h := by
    have h1 : Φ.map ((unTrPre P Fc).Base g) ((unTrPre P Fc).Div φ)
          + (((unTrPre P Fc).degFr φ : ℕ+) : ℕ) • (unTrPre P Fc).Div g
        = Φ.map ((unTrPre P Fc).Base h) ((unTrPre P Fc).Div φ)
          + (((unTrPre P Fc).degFr φ : ℕ+) : ℕ) • (unTrPre P Fc).Div h :=
      congrArg ElemFrobCat.Hom.div e
    rw [hbase, show (unTrPre P Fc).degFr φ = 1 from hφ.1] at h1
    have h2 := add_left_cancel h1
    simpa using h2
  exact (unTrToElem P).map_injective (ElemFrobCat.Hom.ext hbase hdiv hdeg)

/-! ## ★3. ★★★**co-angular は `𝒞^un-tr` では自動**

★★**在庫を引いた結果**(2026-08-17): `prop_1_4_i`(`Prop14.lean`)は
「**始域から出る先がすべて isotropic なら、その射は co-angular**」と言い、
★**`FrobenioidCore` を要求しない**。★`𝒞^un-tr` は全対象が isotropic なので、
★★**すべての射が co-angular** である。

★これで co-angular 系の 5 フィールドと、`Frobenius 型 = 等長な底同型`
(`prop_1_4_i_frobeniusType`)が一気に片付く。 -/

include P in
/-- ★★★**`𝒞^un-tr` ではすべての射が co-angular**。 -/
theorem unTr_coAngular (Fc : FrobenioidCore P) {A B : UnTr P} (f : A ⟶ B) :
    IsCoAngular (unTrPre P Fc) f :=
  prop_1_4_i (unTrPre P Fc) f (fun X _ => unTr_isotropic P Fc X)

include P in
/-- ★★**(iii)(a)** co-angular は合成で閉じる —— ★全射が co-angular なので自明。 -/
theorem unTr_coAngularComp (Fc : FrobenioidCore P) :
    ∀ {A B E : UnTr P} (ψ : A ⟶ B) (φ : B ⟶ E), IsCoAngular (unTrPre P Fc) ψ →
      IsCoAngular (unTrPre P Fc) φ → IsCoAngular (unTrPre P Fc) (ψ ≫ φ) :=
  fun ψ φ _ _ => unTr_coAngular P Fc (ψ ≫ φ)

include P in
/-- ★★**(iii)(b)** —— ★同上。 -/
theorem unTr_coAngularOfPreStep (Fc : FrobenioidCore P) :
    ∀ {A' A : UnTr P} (α : A' ⟶ A), IsCoAngular (unTrPre P Fc) α →
      IsPreStep (unTrPre P Fc) α → ∀ φ : A' ⟶ A, IsCoAngular (unTrPre P Fc) φ :=
  fun _ _ _ φ => unTr_coAngular P Fc φ

include P in
/-- ★★**`𝒞^un-tr` の Frobenius 型 = 等長な底同型**(`Proposition 1.4, (i)` の
「In particular」)。★★これで **Frobenius 型が `Iff.rfl` の 2 条件に潰れる**。 -/
theorem unTr_isFrobeniusType_iff (Fc : FrobenioidCore P) {A B : UnTr P} (φ : A ⟶ B) :
    IsFrobeniusType (unTrPre P Fc) φ
      ↔ IsIsometric (unTrPre P Fc) φ ∧ IsBaseIsomorphism (unTrPre P Fc) φ :=
  prop_1_4_i_frobeniusType (unTrPre P Fc) φ (fun X _ => unTr_isotropic P Fc X)

include P in
/-- ★★**(ii)** 各次数の Frobenius 型射が存在する。

★`𝒞^istr` の証人を送るだけ —— ★**Frobenius 型は「等長 ∧ 底同型」に潰れており、
その 2 つは `Iff.rfl` で商を渡る**。 -/
theorem unTr_frobDegSurj (Fc : FrobenioidCore P) :
    ∀ (A : UnTr P) (n : ℕ+), ∃ (B : UnTr P) (φ : A ⟶ B),
      IsFrobeniusType (unTrPre P Fc) φ ∧ (unTrPre P Fc).degFr φ = n := by
  intro A n
  obtain ⟨B, φ, hft, hdeg⟩ := (istr_frobenioidCore P Fc).frobDegSurj (show Istr P from A) n
  exact ⟨show UnTr P from B, (istrToUnTr P).map φ,
    (unTr_isFrobeniusType_iff P Fc _).mpr ⟨hft.1.2, hft.2⟩, hdeg⟩

/-! ## ★4. pre-step の分解 4 フィールド —— ★★**等長 pre-step が同型**なので潰れる

★★`𝒞^un-tr` は isotropic 型なので、★**等長 pre-step はすべて同型**である。
★したがって `Definition 1.3, (v), (b)(c)` の「co-angular ∘ 等長」分解は
**片方を恒等射に取れば済み**、一意性も**同型の合成で明示的に作れる**。 -/

include P in
/-- ★★**(v)(b)** pre-step は「等長 ∘ co-angular」に分解する —— ★**後半を `𝟙` に取る**。 -/
theorem unTr_preStepFactor (Fc : FrobenioidCore P) :
    ∀ {A B : UnTr P} (φ : A ⟶ B), IsPreStep (unTrPre P Fc) φ →
      ∃ (X : UnTr P) (β : A ⟶ X) (α : X ⟶ B),
        φ = β ≫ α ∧ IsCoAngular (unTrPre P Fc) β ∧ IsPreStep (unTrPre P Fc) β ∧
          IsIsometric (unTrPre P Fc) α ∧ IsPreStep (unTrPre P Fc) α := by
  intro A B φ hφ
  exact ⟨B, φ, 𝟙 B, (Category.comp_id φ).symm, unTr_coAngular P Fc φ, hφ,
    (unTrPre P Fc).Div_id B, isPreStep_id (unTrPre P Fc) B⟩

include P in
/-- ★★**(v)(c)** pre-step は「co-angular ∘ 等長」にも分解する —— ★**前半を `𝟙` に取る**。 -/
theorem unTr_preStepFactor' (Fc : FrobenioidCore P) :
    ∀ {A B : UnTr P} (φ : A ⟶ B), IsPreStep (unTrPre P Fc) φ →
      ∃ (X : UnTr P) (β : A ⟶ X) (α : X ⟶ B),
        φ = β ≫ α ∧ IsIsometric (unTrPre P Fc) β ∧ IsPreStep (unTrPre P Fc) β ∧
          IsCoAngular (unTrPre P Fc) α ∧ IsPreStep (unTrPre P Fc) α := by
  intro A B φ hφ
  exact ⟨A, 𝟙 A, φ, (Category.id_comp φ).symm, (unTrPre P Fc).Div_id A,
    isPreStep_id (unTrPre P Fc) A, unTr_coAngular P Fc φ, hφ⟩

include P in
/-- ★★**(v)(b) の一意性** —— ★**等長 pre-step `α`・`α'` が同型**なので
`γ := α ≫ α'⁻¹` で明示的に作れる。 -/
theorem unTr_preStepFactorUniq (Fc : FrobenioidCore P) :
    ∀ {A B : UnTr P} (X X' : UnTr P) (β : A ⟶ X) (α : X ⟶ B)
      (β' : A ⟶ X') (α' : X' ⟶ B), β ≫ α = β' ≫ α' →
      IsCoAngular (unTrPre P Fc) β → IsPreStep (unTrPre P Fc) β →
      IsIsometric (unTrPre P Fc) α → IsPreStep (unTrPre P Fc) α →
      IsCoAngular (unTrPre P Fc) β' → IsPreStep (unTrPre P Fc) β' →
      IsIsometric (unTrPre P Fc) α' → IsPreStep (unTrPre P Fc) α' →
      ∃ γ : X ≅ X', α' = γ.inv ≫ α ∧ β' = β ≫ γ.hom := by
  intro A B X X' β α β' α' heq _ _ hαi hαs _ _ hα'i hα's
  haveI : IsIso α := unTr_isotropic P Fc X B α hαi hαs
  haveI : IsIso α' := unTr_isotropic P Fc X' B α' hα'i hα's
  refine ⟨asIso α ≪≫ (asIso α').symm, by simp, ?_⟩
  rw [Iso.trans_hom, Iso.symm_hom, asIso_hom, asIso_inv, ← Category.assoc, heq,
    Category.assoc, IsIso.hom_inv_id, Category.comp_id]

include P in
/-- ★★**(v)(c) の一意性** —— ★**等長 pre-step `β`・`β'` が同型**なので
`γ := β⁻¹ ≫ β'` で明示的に作れる。 -/
theorem unTr_preStepFactorUniq' (Fc : FrobenioidCore P) :
    ∀ {A B : UnTr P} (X X' : UnTr P) (β : A ⟶ X) (α : X ⟶ B)
      (β' : A ⟶ X') (α' : X' ⟶ B), β ≫ α = β' ≫ α' →
      IsIsometric (unTrPre P Fc) β → IsPreStep (unTrPre P Fc) β →
      IsCoAngular (unTrPre P Fc) α → IsPreStep (unTrPre P Fc) α →
      IsIsometric (unTrPre P Fc) β' → IsPreStep (unTrPre P Fc) β' →
      IsCoAngular (unTrPre P Fc) α' → IsPreStep (unTrPre P Fc) α' →
      ∃ γ : X ≅ X', α' = γ.inv ≫ α ∧ β' = β ≫ γ.hom := by
  intro A B X X' β α β' α' heq hβi hβs _ _ hβ'i hβ's _ _
  haveI : IsIso β := unTr_isotropic P Fc A X β hβi hβs
  haveI : IsIso β' := unTr_isotropic P Fc A X' β' hβ'i hβ's
  refine ⟨(asIso β).symm ≪≫ asIso β', ?_, by simp⟩
  rw [Iso.trans_inv, Iso.symm_inv, asIso_hom, asIso_inv, Category.assoc, heq,
    ← Category.assoc, IsIso.inv_hom_id, Category.id_comp]

/-! ## ★5. `𝒪^▷` の全単射(3 フィールド)

★★`otriBase` は ★**忠実性で直接出る** —— 条件は `𝔽_Φ` の中の等式に翻訳され、
`Base φ = Base φ'` がそのまま効く。 -/

include P in
/-- ★★**(iii)(c)** 全単射は `Base(φ)` にしか依らない。

★★`φ ≫ β = α ≫ φ` を `𝔽_Φ` に落とすと **`Φ.map (Base φ) (Div β) = Div α`** になる。
★これは `Base φ` にしか依らないので、`φ'` でも同じ式が成り立つ。 -/
theorem unTr_otriBase (Fc : FrobenioidCore P) :
    ∀ {A B : UnTr P} (φ φ' : A ⟶ B), IsCoAngular (unTrPre P Fc) φ →
      IsPreStep (unTrPre P Fc) φ → IsCoAngular (unTrPre P Fc) φ' →
      IsPreStep (unTrPre P Fc) φ' →
      (unTrPre P Fc).Base φ = (unTrPre P Fc).Base φ' →
      ∀ α ∈ OTri (unTrPre P Fc) A, ∀ β ∈ OTri (unTrPre P Fc) B,
        (φ ≫ (β : B ⟶ B) : A ⟶ B) = (α : A ⟶ A) ≫ φ →
        (φ' ≫ (β : B ⟶ B) : A ⟶ B) = (α : A ⟶ A) ≫ φ' := by
  intro A B φ φ' _ hs _ hs' hbase α hα β hβ h
  letI := isCancelAdd_of_isIntegralMonoid (Φ.val ((unTrPre P Fc).toElem.obj A).base)
    (P.divisorial _).1.1
  have hαb : (unTrPre P Fc).Base (α : A ⟶ A) = 𝟙 _ := by
    have h0 : (unTrPre P Fc).Base (α : A ⟶ A) = (unTrPre P Fc).Base (𝟙 A) := hα.1
    rw [h0, (unTrPre P Fc).Base_id]
  have hβb : (unTrPre P Fc).Base (β : B ⟶ B) = 𝟙 _ := by
    have h0 : (unTrPre P Fc).Base (β : B ⟶ B) = (unTrPre P Fc).Base (𝟙 B) := hβ.1
    rw [h0, (unTrPre P Fc).Base_id]
  -- ★核心の式を取り出す
  have key : Φ.map ((unTrPre P Fc).Base φ) ((unTrPre P Fc).Div (β : B ⟶ B))
      = (unTrPre P Fc).Div (α : A ⟶ A) := by
    have hE : (unTrToElem P).map φ ≫ (unTrToElem P).map (β : B ⟶ B)
        = (unTrToElem P).map (α : A ⟶ A) ≫ (unTrToElem P).map φ := by
      rw [← (unTrToElem P).map_comp, ← (unTrToElem P).map_comp, h]
    have hd : Φ.map ((unTrPre P Fc).Base φ) ((unTrPre P Fc).Div (β : B ⟶ B))
          + (((unTrPre P Fc).degFr (β : B ⟶ B) : ℕ+) : ℕ) • (unTrPre P Fc).Div φ
        = Φ.map ((unTrPre P Fc).Base (α : A ⟶ A)) ((unTrPre P Fc).Div φ)
          + (((unTrPre P Fc).degFr φ : ℕ+) : ℕ) • (unTrPre P Fc).Div (α : A ⟶ A) :=
      congrArg ElemFrobCat.Hom.div hE
    rw [show (unTrPre P Fc).degFr (β : B ⟶ B) = 1 from hβ.2,
      show (unTrPre P Fc).degFr φ = 1 from hs.1, hαb, Φ.map_id] at hd
    simp only [PNat.one_coe, one_smul] at hd
    rw [add_comm ((unTrPre P Fc).Div φ) ((unTrPre P Fc).Div (α : A ⟶ A))] at hd
    exact add_right_cancel hd
  refine (unTrToElem P).map_injective ?_
  rw [(unTrToElem P).map_comp, (unTrToElem P).map_comp]
  refine ElemFrobCat.Hom.ext ?_ ?_ ?_
  · show (unTrPre P Fc).Base φ' ≫ (unTrPre P Fc).Base (β : B ⟶ B)
      = (unTrPre P Fc).Base (α : A ⟶ A) ≫ (unTrPre P Fc).Base φ'
    rw [hαb, hβb, Category.comp_id, Category.id_comp]
  · show Φ.map ((unTrPre P Fc).Base φ') ((unTrPre P Fc).Div (β : B ⟶ B))
        + (((unTrPre P Fc).degFr (β : B ⟶ B) : ℕ+) : ℕ) • (unTrPre P Fc).Div φ'
      = Φ.map ((unTrPre P Fc).Base (α : A ⟶ A)) ((unTrPre P Fc).Div φ')
        + (((unTrPre P Fc).degFr φ' : ℕ+) : ℕ) • (unTrPre P Fc).Div (α : A ⟶ A)
    rw [show (unTrPre P Fc).degFr (β : B ⟶ B) = 1 from hβ.2,
      show (unTrPre P Fc).degFr φ' = 1 from hs'.1, hαb, Φ.map_id, ← hbase, key]
    simp only [PNat.one_coe, one_smul]
    exact add_comm _ _
  · show (unTrPre P Fc).degFr (β : B ⟶ B) * (unTrPre P Fc).degFr φ'
      = (unTrPre P Fc).degFr φ' * (unTrPre P Fc).degFr (α : A ⟶ A)
    rw [show (unTrPre P Fc).degFr (β : B ⟶ B) = 1 from hβ.2,
      show (unTrPre P Fc).degFr φ' = 1 from hs'.1,
      show (unTrPre P Fc).degFr (α : A ⟶ A) = 1 from hα.2]

include P in
/-- ★★**(iii)(c)** co-angular pre-step は `𝒪^▷` の全単射を誘導する(順向き)。

★★**存在は `𝒞^istr` から送る**(`istrToUnTr` は full なので持ち上げられる)。
★★**一意性は `𝒞^un-tr` が totally epimorphic** であることから ——
`φ ≫ β₁ = φ ≫ β₂` を `φ` が epi なので消せる。 -/
theorem unTr_otriFwd (Fc : FrobenioidCore P) :
    ∀ {A B : UnTr P} (φ : A ⟶ B), IsCoAngular (unTrPre P Fc) φ →
      IsPreStep (unTrPre P Fc) φ →
      ∀ α ∈ OTri (unTrPre P Fc) A, ∃! β : End B, β ∈ OTri (unTrPre P Fc) B ∧
        (φ ≫ (β : B ⟶ B) : A ⟶ B) = (α : A ⟶ A) ≫ φ := by
  intro A B φ _ hs α hα
  obtain ⟨φ₀, rfl⟩ := (istrToUnTr P).map_surjective φ
  obtain ⟨α₀, rfl⟩ := (istrToUnTr P).map_surjective (α : A ⟶ A)
  obtain ⟨β₀, ⟨hβ₀mem, hβ₀eq⟩, -⟩ :=
    (istr_frobenioidCore P Fc).otriFwd φ₀ (istr_coAngular P Fc φ₀) hs α₀ hα
  haveI : Epi ((istrToUnTr P).map φ₀) := unTr_totEpi P _ _ _
  refine ⟨(istrToUnTr P).map β₀, ⟨hβ₀mem, ?_⟩, ?_⟩
  · exact congrArg (fun t => (istrToUnTr P).map t) hβ₀eq
  · rintro β ⟨-, hβeq⟩
    refine (cancel_epi ((istrToUnTr P).map φ₀)).mp ?_
    exact hβeq.trans (congrArg (fun t => (istrToUnTr P).map t) hβ₀eq).symm

include P in
/-- ★★**(iii)(c)** 逆向き。★手は `unTr_otriFwd` と同じ。 -/
theorem unTr_otriBwd (Fc : FrobenioidCore P) :
    ∀ {A B : UnTr P} (φ : A ⟶ B), IsCoAngular (unTrPre P Fc) φ →
      IsPreStep (unTrPre P Fc) φ →
      ∀ β ∈ OTri (unTrPre P Fc) B, ∃! α : End A, α ∈ OTri (unTrPre P Fc) A ∧
        (φ ≫ (β : B ⟶ B) : A ⟶ B) = (α : A ⟶ A) ≫ φ := by
  intro A B φ _ hs β hβ
  obtain ⟨φ₀, rfl⟩ := (istrToUnTr P).map_surjective φ
  obtain ⟨β₀, rfl⟩ := (istrToUnTr P).map_surjective (β : B ⟶ B)
  obtain ⟨α₀, ⟨hα₀mem, hα₀eq⟩, -⟩ :=
    (istr_frobenioidCore P Fc).otriBwd φ₀ (istr_coAngular P Fc φ₀) hs β₀ hβ
  haveI : Mono ((istrToUnTr P).map φ₀) := unTr_preStepMono P Fc _ hs
  refine ⟨(istrToUnTr P).map α₀, ⟨hα₀mem, ?_⟩, ?_⟩
  · exact congrArg (fun t => (istrToUnTr P).map t) hα₀eq
  · rintro α ⟨-, hαeq⟩
    refine (cancel_mono ((istrToUnTr P).map φ₀)).mp ?_
    exact hαeq.symm.trans (congrArg (fun t => (istrToUnTr P).map t) hα₀eq)

/-! ## ★6. 持ち上げて送る 2 フィールド -/

include P in
/-- ★★**(ii)** Frobenius 型射の本質的一意性。

★`istrToUnTr` は **full** なので `𝒞^istr` へ持ち上げ、`istr_frobDegUniq` を当てて送り返す。
★★**Frobenius 型が持ち上がるのは、co-angular が `𝒞^istr` でも自動**だから
(`istr_coAngular`)——等長と底同型は `Iff.rfl` で渡る。 -/
theorem unTr_frobDegUniq (Fc : FrobenioidCore P) :
    ∀ (A B E : UnTr P) (φ : A ⟶ B) (ψ : A ⟶ E), IsFrobeniusType (unTrPre P Fc) φ →
      IsFrobeniusType (unTrPre P Fc) ψ →
      (unTrPre P Fc).degFr φ = (unTrPre P Fc).degFr ψ →
      ∃ β : B ⟶ E, IsIso β ∧ φ ≫ β = ψ := by
  intro A B E φ ψ hφ hψ hdeg
  obtain ⟨φ₀, rfl⟩ := (istrToUnTr P).map_surjective φ
  obtain ⟨ψ₀, rfl⟩ := (istrToUnTr P).map_surjective ψ
  obtain ⟨β₀, hβiso, hβeq⟩ := (istr_frobenioidCore P Fc).frobDegUniq
    (show Istr P from A) (show Istr P from B) (show Istr P from E) φ₀ ψ₀
    ⟨⟨istr_coAngular P Fc φ₀, hφ.1.2⟩, hφ.2⟩
    ⟨⟨istr_coAngular P Fc ψ₀, hψ.1.2⟩, hψ.2⟩ hdeg
  haveI := hβiso
  refine ⟨(istrToUnTr P).map β₀, ?_, congrArg (fun t => (istrToUnTr P).map t) hβeq⟩
  exact ⟨⟨((istrToUnTr P).mapIso (asIso β₀)).inv,
    ((istrToUnTr P).mapIso (asIso β₀)).hom_inv_id,
    ((istrToUnTr P).mapIso (asIso β₀)).inv_hom_id⟩⟩

include P in
/-- ★★**(i)(a)** 底の圏への全射性。

★`𝒞^istr` の証人を送る。★**`Frobenius-trivial` の族 `ζ` は
`End` の単系準同型として押し出せる**(`istrToUnTr` が関手だから)。 -/
theorem unTr_baseSurj (Fc : FrobenioidCore P) :
    ∀ Y : D, ∃ A : UnTr P, IsFrobeniusTrivial (unTrPre P Fc) A ∧
      Nonempty (((unTrPre P Fc).toElem.obj A).base ≅ Y) := by
  intro Y
  obtain ⟨A, ⟨ζ, hdeg, hprop⟩, he⟩ := (istr_frobenioidCore P Fc).baseSurj Y
  refine ⟨show UnTr P from A, ⟨?_, ?_⟩, he⟩
  · exact
      { toFun := fun n => (istrToUnTr P).map (ζ n)
        map_one' := by
          show (istrToUnTr P).map (ζ 1) = 𝟙 _
          rw [map_one]
          exact (istrToUnTr P).map_id _
        map_mul' := fun x y => by
          show (istrToUnTr P).map (ζ (x * y)) = _
          rw [map_mul]
          exact ((istrToUnTr P).map_comp (ζ y) (ζ x)).symm }
  · refine ⟨fun n => hdeg n, fun n => ⟨(hprop n).1, ?_⟩⟩
    exact (unTr_isFrobeniusType_iff P Fc _).mpr ⟨(hprop n).2.1.2, (hprop n).2.2⟩

/-! ## ★7. pull-back 射は `𝒞^istr` から `𝒞^un-tr` へ渡る

★★**単射性は `𝒞^un-tr` では自動**である —— 忠実性より射は `(Base, Div, degFr)` で
決まり、`Base` は仮定から、`degFr` は `ℕ+` の簡約から、`Div` は
★**pull-back 射が linear(`degFr = 1`)**なので `1 • (−)` の簡約から一致する。

★★**「pull-back は linear」を先に使うのが要点**である —— これがないと
`n • Div f₁ = n • Div f₂` から `Div f₁ = Div f₂` を出すために
`Gp Φ` の捻れ無しが要る(それも真だが、遠回りである)。

★★全射性は `𝒞^istr` の証人を送るだけ。 -/

include P in
/-- ★★★**pull-back 射は商へ渡る**。 -/
theorem unTr_isPullBack_of_istr (Fc : FrobenioidCore P) {A B : Istr P} (φ₀ : A ⟶ B)
    (hpb : IsPullBack (istrPre P Fc) φ₀) :
    IsPullBack (unTrPre P Fc) ((istrToUnTr P).map φ₀) := by
  have hlin : (istrPre P Fc).degFr φ₀ = 1 := ((istr_frobenioidCore P Fc).pullBackLB φ₀ hpb).2
  intro X
  constructor
  · -- ★単射性
    intro f₁ f₂ hf
    have hf' := Subtype.ext_iff.mp hf
    have hcomp : f₁ ≫ (istrToUnTr P).map φ₀ = f₂ ≫ (istrToUnTr P).map φ₀ :=
      congrArg Prod.fst hf'
    have hbase : (unTrPre P Fc).Base f₁ = (unTrPre P Fc).Base f₂ := congrArg Prod.snd hf'
    letI : IsCancelAdd (Φ.val ((unTrToElem P).obj X).base) :=
      isCancelAdd_of_isIntegralMonoid _ (P.divisorial ((unTrToElem P).obj X).base).1.1
    have hE : (unTrToElem P).map f₁ ≫ (unTrToElem P).map ((istrToUnTr P).map φ₀)
        = (unTrToElem P).map f₂ ≫ (unTrToElem P).map ((istrToUnTr P).map φ₀) := by
      rw [← (unTrToElem P).map_comp, ← (unTrToElem P).map_comp, hcomp]
    refine (unTrToElem P).map_injective (ElemFrobCat.Hom.ext hbase ?_ ?_)
    · have hd : Φ.map ((unTrPre P Fc).Base f₁)
            ((unTrPre P Fc).Div ((istrToUnTr P).map φ₀))
            + (((istrPre P Fc).degFr φ₀ : ℕ+) : ℕ) • (unTrPre P Fc).Div f₁
          = Φ.map ((unTrPre P Fc).Base f₂)
            ((unTrPre P Fc).Div ((istrToUnTr P).map φ₀))
            + (((istrPre P Fc).degFr φ₀ : ℕ+) : ℕ) • (unTrPre P Fc).Div f₂ :=
        congrArg ElemFrobCat.Hom.div hE
      rw [hbase, hlin, PNat.one_coe, one_smul, one_smul] at hd
      -- ★instance 合成ではなく**明示適用**で通す(defeq 検査になるので食い違いを跨げる)
      exact @add_left_cancel _ _ (isCancelAdd_of_isIntegralMonoid
        (Φ.val ((unTrToElem P).obj X).base) (P.divisorial _).1.1).toIsLeftCancelAdd _ _ _ hd
    · have hg : (istrPre P Fc).degFr φ₀ * (unTrPre P Fc).degFr f₁
          = (istrPre P Fc).degFr φ₀ * (unTrPre P Fc).degFr f₂ :=
        congrArg ElemFrobCat.Hom.deg hE
      exact mul_left_cancel hg
  · -- ★全射性
    rintro ⟨⟨g, hh⟩, hgh⟩
    obtain ⟨g₀, rfl⟩ := (istrToUnTr P).map_surjective g
    obtain ⟨f₀, hf₀⟩ := (hpb (show Istr P from X)).2 ⟨(g₀, hh), hgh⟩
    have hpair : (f₀ ≫ φ₀, (istrPre P Fc).Base f₀) = (g₀, hh) := Subtype.ext_iff.mp hf₀
    have h1 : (istrToUnTr P).map f₀ ≫ (istrToUnTr P).map φ₀ = (istrToUnTr P).map g₀ :=
      congrArg (fun t => (istrToUnTr P).map t) (congrArg Prod.fst hpair)
    have h2 : (unTrPre P Fc).Base ((istrToUnTr P).map f₀) = hh := congrArg Prod.snd hpair
    exact ⟨(istrToUnTr P).map f₀, Subtype.ext (congrArg₂ Prod.mk h1 h2)⟩

include P in
/-- ★★**(iv)(a)** 任意射の 3 分解 —— `𝒞^istr` の分解を送るだけ。

★★**3 つの型がすべて渡る**ことが要る: Frobenius 型(★co-angular が自動なので
「等長 ∧ 底同型」に潰れる)、pre-step(`Iff.rfl`)、pull-back(★上の `unTr_isPullBack_of_istr`)。 -/
theorem unTr_arbFactor (Fc : FrobenioidCore P) :
    ∀ {A B : UnTr P} (φ : A ⟶ B), ∃ (X Y : UnTr P) (γ : A ⟶ X) (β : X ⟶ Y) (α : Y ⟶ B),
      φ = γ ≫ β ≫ α ∧ IsFrobeniusType (unTrPre P Fc) γ ∧ IsPreStep (unTrPre P Fc) β ∧
        IsPullBack (unTrPre P Fc) α := by
  intro A B φ
  obtain ⟨φ₀, rfl⟩ := (istrToUnTr P).map_surjective φ
  obtain ⟨X, Y, γ, β, α, heq, hγ, hβ, hα⟩ := (istr_frobenioidCore P Fc).arbFactor φ₀
  exact ⟨show UnTr P from X, show UnTr P from Y, (istrToUnTr P).map γ,
    (istrToUnTr P).map β, (istrToUnTr P).map α,
    congrArg (fun t => (istrToUnTr P).map t) heq,
    (unTr_isFrobeniusType_iff P Fc _).mpr ⟨hγ.1.2, hγ.2⟩, hβ,
    unTr_isPullBack_of_istr P Fc α hα⟩

/-! ## ★7-b. ★★★**pull-back 射の逆向きの移送**(2026-08-17)

★★**どの持ち上げでも pull-back になる。**★筋は「**余分な部分が同型に潰れる**」:

1. `𝒞^istr` の `arbFactor` で `α₀ = γ ≫ β ≫ α'`
2. ★★**`Base (γ ≫ β)` は同型**(Frobenius 型も pre-step も底同型)——★ここが鍵
3. `IsPullBack f` の**全射性**を対 `(map α', inv (Base (γ≫β)))` に当てて `d` を得る
4. 両側の**単射性**で `map(γ≫β)` が `𝒞^un-tr` の同型と分かる
5. ★**同型 ⟹ 等長 pre-step ⟹ `𝒞^istr` は isotropic 型なので `γ ≫ β` も同型**
6. `IsPullBack.comp` ＋ `isPullBack_of_isIso`(★どちらも在庫)で終わり

★★**一度「unit-equivalence で不変でないから駄目」と判断したが誤りだった** ——
★**不変性は要らず、分解して潰せばよい**。 -/

include P in
/-- ★★★**`𝒞^un-tr` の pull-back 射の持ち上げは `𝒞^istr` の pull-back 射**。 -/
theorem unTr_isPullBack_to (Fc : FrobenioidCore P) {A B : Istr P} (α₀ : A ⟶ B)
    (h : IsPullBack (unTrPre P Fc) ((istrToUnTr P).map α₀)) :
    IsPullBack (istrPre P Fc) α₀ := by
  obtain ⟨Y₀, Y, γ, β, α', heq, hγ, hβ, hα'⟩ := (istr_frobenioidCore P Fc).arbFactor α₀
  haveI hbγ : IsIso ((istrPre P Fc).Base γ) := hγ.2
  haveI hbβ : IsIso ((istrPre P Fc).Base β) := hβ.2
  haveI hbc : IsIso ((istrPre P Fc).Base (γ ≫ β)) := by
    rw [(istrPre P Fc).Base_comp]
    exact IsIso.comp_isIso' hbγ hbβ
  have heq' : α₀ = (γ ≫ β) ≫ α' := heq.trans (Category.assoc _ _ _).symm
  have hfc : (istrToUnTr P).map α₀
      = (istrToUnTr P).map (γ ≫ β) ≫ (istrToUnTr P).map α' :=
    congrArg (fun t => (istrToUnTr P).map t) heq'
  -- ★★底の同型を**不透明な局所変数**にする(`set` だと instance 探索が通らない)
  obtain ⟨u, hu⟩ : ∃ u : ((unTrPre P Fc).toElem.obj (show UnTr P from A)).base
      ⟶ ((unTrPre P Fc).toElem.obj (show UnTr P from Y)).base,
      u = (unTrPre P Fc).Base ((istrToUnTr P).map (γ ≫ β)) := ⟨_, rfl⟩
  haveI huiso : IsIso u := by rw [hu]; exact hbc
  have hpb' : IsPullBack (unTrPre P Fc) ((istrToUnTr P).map α') :=
    unTr_isPullBack_of_istr P Fc α' hα'
  -- ★手 3: `d` を作る
  have hbα₀ : (unTrPre P Fc).Base ((istrToUnTr P).map α₀)
      = u ≫ (unTrPre P Fc).Base ((istrToUnTr P).map α') := by
    rw [hu, hfc]
    exact (unTrPre P Fc).Base_comp _ _
  have hba' : (unTrPre P Fc).Base ((istrToUnTr P).map α')
      = inv u ≫ (unTrPre P Fc).Base ((istrToUnTr P).map α₀) :=
    calc (unTrPre P Fc).Base ((istrToUnTr P).map α')
        = 𝟙 _ ≫ (unTrPre P Fc).Base ((istrToUnTr P).map α') := (Category.id_comp _).symm
      _ = (inv u ≫ u) ≫ (unTrPre P Fc).Base ((istrToUnTr P).map α') := by
          rw [IsIso.inv_hom_id]
          rfl
      _ = inv u ≫ u ≫ (unTrPre P Fc).Base ((istrToUnTr P).map α') := Category.assoc _ _ _
      _ = inv u ≫ (unTrPre P Fc).Base ((istrToUnTr P).map α₀) :=
          congrArg (fun t => inv u ≫ t) hbα₀.symm
  obtain ⟨d, hd⟩ := (h (show UnTr P from Y)).2
    ⟨((istrToUnTr P).map α', @inv _ _ _ _ u huiso), hba'⟩
  have hdpair : (d ≫ (istrToUnTr P).map α₀, (unTrPre P Fc).Base d)
      = ((istrToUnTr P).map α', @inv _ _ _ _ u huiso) := Subtype.ext_iff.mp hd
  have hd1 : d ≫ (istrToUnTr P).map α₀ = (istrToUnTr P).map α' := congrArg Prod.fst hdpair
  have hd2 : (unTrPre P Fc).Base d = @inv _ _ _ _ u huiso := congrArg Prod.snd hdpair
  -- ★手 4: `map (γ ≫ β)` が同型
  have hcd : (istrToUnTr P).map (γ ≫ β) ≫ d = 𝟙 _ := by
    refine (h (show UnTr P from A)).1 (Subtype.ext (congrArg₂ Prod.mk ?_ ?_))
    · show ((istrToUnTr P).map (γ ≫ β) ≫ d) ≫ (istrToUnTr P).map α₀
        = 𝟙 _ ≫ (istrToUnTr P).map α₀
      calc ((istrToUnTr P).map (γ ≫ β) ≫ d) ≫ (istrToUnTr P).map α₀
          = (istrToUnTr P).map (γ ≫ β) ≫ d ≫ (istrToUnTr P).map α₀ := Category.assoc _ _ _
        _ = (istrToUnTr P).map (γ ≫ β) ≫ (istrToUnTr P).map α' :=
            congrArg (fun t => (istrToUnTr P).map (γ ≫ β) ≫ t) hd1
        _ = (istrToUnTr P).map α₀ := hfc.symm
        _ = 𝟙 _ ≫ (istrToUnTr P).map α₀ := (Category.id_comp _).symm
    · show (unTrPre P Fc).Base ((istrToUnTr P).map (γ ≫ β) ≫ d)
        = (unTrPre P Fc).Base (𝟙 _)
      calc (unTrPre P Fc).Base ((istrToUnTr P).map (γ ≫ β) ≫ d)
          = (unTrPre P Fc).Base ((istrToUnTr P).map (γ ≫ β)) ≫ (unTrPre P Fc).Base d :=
            (unTrPre P Fc).Base_comp _ _
        _ = u ≫ inv u := congrArg₂ (fun x y => x ≫ y) hu.symm hd2
        _ = 𝟙 _ := IsIso.hom_inv_id u
        _ = (unTrPre P Fc).Base (𝟙 _) := ((unTrPre P Fc).Base_id _).symm
  have hdc : d ≫ (istrToUnTr P).map (γ ≫ β) = 𝟙 _ := by
    refine (hpb' (show UnTr P from Y)).1 (Subtype.ext (congrArg₂ Prod.mk ?_ ?_))
    · show (d ≫ (istrToUnTr P).map (γ ≫ β)) ≫ (istrToUnTr P).map α'
        = 𝟙 _ ≫ (istrToUnTr P).map α'
      calc (d ≫ (istrToUnTr P).map (γ ≫ β)) ≫ (istrToUnTr P).map α'
          = d ≫ (istrToUnTr P).map (γ ≫ β) ≫ (istrToUnTr P).map α' := Category.assoc _ _ _
        _ = d ≫ (istrToUnTr P).map α₀ :=
            congrArg (fun t => d ≫ t) hfc.symm
        _ = (istrToUnTr P).map α' := hd1
        _ = 𝟙 _ ≫ (istrToUnTr P).map α' := (Category.id_comp _).symm
    · show (unTrPre P Fc).Base (d ≫ (istrToUnTr P).map (γ ≫ β))
        = (unTrPre P Fc).Base (𝟙 _)
      calc (unTrPre P Fc).Base (d ≫ (istrToUnTr P).map (γ ≫ β))
          = (unTrPre P Fc).Base d ≫ (unTrPre P Fc).Base ((istrToUnTr P).map (γ ≫ β)) :=
            (unTrPre P Fc).Base_comp _ _
        _ = inv u ≫ u := congrArg₂ (fun x y => x ≫ y) hd2 hu.symm
        _ = 𝟙 _ := IsIso.inv_hom_id u
        _ = (unTrPre P Fc).Base (𝟙 _) := ((unTrPre P Fc).Base_id _).symm
  haveI hciso : IsIso ((istrToUnTr P).map (γ ≫ β)) := ⟨d, hcd, hdc⟩
  -- ★手 5: `γ ≫ β` は `𝒞^istr` の同型
  haveI hE : IsIso ((unTrToElem P).map ((istrToUnTr P).map (γ ≫ β))) :=
    ⟨⟨((unTrToElem P).mapIso (asIso ((istrToUnTr P).map (γ ≫ β)))).inv,
      ((unTrToElem P).mapIso (asIso ((istrToUnTr P).map (γ ≫ β)))).hom_inv_id,
      ((unTrToElem P).mapIso (asIso ((istrToUnTr P).map (γ ≫ β)))).inv_hom_id⟩⟩
  obtain ⟨-, hu0, hd0⟩ := (ElemFrobCat.isIso_iff
    ((unTrToElem P).map ((istrToUnTr P).map (γ ≫ β)))).mp hE
  have hdiv0 : (istrPre P Fc).Div (γ ≫ β) = 0 :=
    (P.divisorial ((istrPre P Fc).toElem.obj A).base).2 _ hu0
  haveI : IsIso (γ ≫ β) :=
    istr_isotropic P Fc A (show Istr P from Y) (γ ≫ β) hdiv0 ⟨hd0, hbc⟩
  -- ★手 6
  rw [heq']
  exact IsPullBack.comp (istrPre P Fc) (isPullBack_of_isIso (istrPre P Fc) (γ ≫ β)) hα'

include P in
/-- ★★**(iv)(b)** pull-back 射は LB-invertible かつ linear。

★★上の逆向き移送で `𝒞^istr` へ降ろし、`istr_pullBackLB` を当てて戻すだけ。
★co-angular は `𝒞^un-tr` では自動なので、`hlb.2`(等長)だけ拾えばよい。 -/
theorem unTr_pullBackLB (Fc : FrobenioidCore P) :
    ∀ {A B : UnTr P} (α : A ⟶ B), IsPullBack (unTrPre P Fc) α →
      IsLBInvertible (unTrPre P Fc) α ∧ IsLinear (unTrPre P Fc) α := by
  intro A B α hα
  obtain ⟨α₀, rfl⟩ := (istrToUnTr P).map_surjective α
  obtain ⟨hlb, hlin⟩ :=
    (istr_frobenioidCore P Fc).pullBackLB α₀ (unTr_isPullBack_to P Fc α₀ hα)
  exact ⟨⟨unTr_coAngular P Fc _, hlb.2⟩, hlin⟩

include P in
/-- ★★**(iv)(a) の一意性** —— ★**持ち上げず、`𝒞^un-tr` の中で直接**組む。

★★**筋**(2026-08-17):
1. ★`pullBackLB` で `α`・`α'` が linear ⟹ **次数は `γ` が全部担う**
2. `frobDegUniq` で同型 `ε : X ≅ X'`、`γ ≫ ε = γ'`
3. ★**`γ` は epi** なので消して `β ≫ α = (ε ≫ β') ≫ α'`
4. ★`IsPullBack α'` / `IsPullBack α` の**全射性**で `δ`・`δ⁻¹` を作り、
   **単射性**で互いに逆と分かる
5. 残る式も `IsPullBack α'` の**単射性**で出る

★★**依存の順序が `pullBackLB → arbFactorUniq`** である(1 でそれが要る)。 -/
theorem unTr_arbFactorUniq (Fc : FrobenioidCore P) :
    ∀ {A B : UnTr P} (X Y X' Y' : UnTr P)
      (γ : A ⟶ X) (β : X ⟶ Y) (α : Y ⟶ B) (γ' : A ⟶ X') (β' : X' ⟶ Y') (α' : Y' ⟶ B),
      γ ≫ β ≫ α = γ' ≫ β' ≫ α' →
      IsFrobeniusType (unTrPre P Fc) γ → IsPreStep (unTrPre P Fc) β →
      IsPullBack (unTrPre P Fc) α →
      IsFrobeniusType (unTrPre P Fc) γ' → IsPreStep (unTrPre P Fc) β' →
      IsPullBack (unTrPre P Fc) α' →
      ∃ (δ : Y ≅ Y') (ε : X ≅ X'),
        α' = δ.inv ≫ α ∧ β' = ε.inv ≫ β ≫ δ.hom ∧ γ' = γ ≫ ε.hom := by
  intro A B X Y X' Y' γ β α γ' β' α' heq hγ hβ hα hγ' hβ' hα'
  have hdα : (unTrPre P Fc).degFr α = 1 := (unTr_pullBackLB P Fc α hα).2
  have hdα' : (unTrPre P Fc).degFr α' = 1 := (unTr_pullBackLB P Fc α' hα').2
  -- ★手 1: 次数は `γ` が担う
  have hdγ : (unTrPre P Fc).degFr γ = (unTrPre P Fc).degFr γ' := by
    have e1 : (unTrPre P Fc).degFr (γ ≫ β ≫ α) = (unTrPre P Fc).degFr γ := by
      rw [(unTrPre P Fc).degFr_comp, (unTrPre P Fc).degFr_comp, hdα,
        show (unTrPre P Fc).degFr β = 1 from hβ.1, one_mul, one_mul]
    have e2 : (unTrPre P Fc).degFr (γ' ≫ β' ≫ α') = (unTrPre P Fc).degFr γ' := by
      rw [(unTrPre P Fc).degFr_comp, (unTrPre P Fc).degFr_comp, hdα',
        show (unTrPre P Fc).degFr β' = 1 from hβ'.1, one_mul, one_mul]
    rw [← e1, ← e2, heq]
  -- ★手 2: `ε`
  obtain ⟨ε₀, hε₀iso, hε₀⟩ := unTr_frobDegUniq P Fc A X X' γ γ' hγ hγ' hdγ
  haveI := hε₀iso
  -- ★手 3: `γ` を消す
  haveI hγepi : Epi γ := unTr_totEpi P _ _ _
  have hβ'' : IsPreStep (unTrPre P Fc) (ε₀ ≫ β') :=
    IsPreStep.comp (unTrPre P Fc) (isPreStep_of_isIso (unTrPre P Fc) ε₀) hβ'
  have hcancel : β ≫ α = (ε₀ ≫ β') ≫ α' := by
    refine (cancel_epi γ).mp ?_
    calc γ ≫ β ≫ α = γ' ≫ β' ≫ α' := heq
      _ = (γ ≫ ε₀) ≫ β' ≫ α' := by rw [hε₀]
      _ = γ ≫ (ε₀ ≫ β') ≫ α' := by simp only [Category.assoc]
  -- ★手 4: `δ` を作る
  haveI hbβ : IsIso ((unTrPre P Fc).Base β) := hβ.2
  haveI hbβ'' : IsIso ((unTrPre P Fc).Base (ε₀ ≫ β')) := hβ''.2
  obtain ⟨hh, hhdef⟩ : ∃ hh : ((unTrPre P Fc).toElem.obj Y).base
      ⟶ ((unTrPre P Fc).toElem.obj Y').base,
      hh = inv ((unTrPre P Fc).Base β) ≫ (unTrPre P Fc).Base (ε₀ ≫ β') := ⟨_, rfl⟩
  haveI hhiso : IsIso hh := by rw [hhdef]; infer_instance
  have hbaseα : (unTrPre P Fc).Base α = hh ≫ (unTrPre P Fc).Base α' := by
    have e : (unTrPre P Fc).Base β ≫ (unTrPre P Fc).Base α
        = (unTrPre P Fc).Base (ε₀ ≫ β') ≫ (unTrPre P Fc).Base α' := by
      rw [← (unTrPre P Fc).Base_comp, ← (unTrPre P Fc).Base_comp, hcancel]
    rw [hhdef, Category.assoc, ← e, ← Category.assoc, IsIso.inv_hom_id, Category.id_comp]
  have hbaseα' : (unTrPre P Fc).Base α'
      = @inv _ _ _ _ hh hhiso ≫ (unTrPre P Fc).Base α := by
    rw [hbaseα, ← Category.assoc, IsIso.inv_hom_id, Category.id_comp]
  obtain ⟨δ₀, hδ₀⟩ := (hα' Y).2 ⟨(α, hh), hbaseα⟩
  have hδ₀p : (δ₀ ≫ α', (unTrPre P Fc).Base δ₀) = (α, hh) := Subtype.ext_iff.mp hδ₀
  have hδ₀1 : δ₀ ≫ α' = α := congrArg Prod.fst hδ₀p
  have hδ₀2 : (unTrPre P Fc).Base δ₀ = hh := congrArg Prod.snd hδ₀p
  obtain ⟨δ₁, hδ₁⟩ := (hα Y').2 ⟨(α', @inv _ _ _ _ hh hhiso), hbaseα'⟩
  have hδ₁p : (δ₁ ≫ α, (unTrPre P Fc).Base δ₁) = (α', @inv _ _ _ _ hh hhiso) :=
    Subtype.ext_iff.mp hδ₁
  have hδ₁1 : δ₁ ≫ α = α' := congrArg Prod.fst hδ₁p
  have hδ₁2 : (unTrPre P Fc).Base δ₁ = @inv _ _ _ _ hh hhiso := congrArg Prod.snd hδ₁p
  -- ★手 4b: 互いに逆
  have hδid : δ₀ ≫ δ₁ = 𝟙 Y := by
    refine (hα Y).1 (Subtype.ext (congrArg₂ Prod.mk ?_ ?_))
    · show (δ₀ ≫ δ₁) ≫ α = 𝟙 Y ≫ α
      calc (δ₀ ≫ δ₁) ≫ α = δ₀ ≫ δ₁ ≫ α := Category.assoc _ _ _
        _ = δ₀ ≫ α' := congrArg (fun t => δ₀ ≫ t) hδ₁1
        _ = α := hδ₀1
        _ = 𝟙 Y ≫ α := (Category.id_comp _).symm
    · show (unTrPre P Fc).Base (δ₀ ≫ δ₁) = (unTrPre P Fc).Base (𝟙 Y)
      calc (unTrPre P Fc).Base (δ₀ ≫ δ₁)
          = (unTrPre P Fc).Base δ₀ ≫ (unTrPre P Fc).Base δ₁ := (unTrPre P Fc).Base_comp _ _
        _ = hh ≫ @inv _ _ _ _ hh hhiso := congrArg₂ (fun x y => x ≫ y) hδ₀2 hδ₁2
        _ = 𝟙 _ := IsIso.hom_inv_id hh
        _ = (unTrPre P Fc).Base (𝟙 Y) := ((unTrPre P Fc).Base_id _).symm
  have hδid' : δ₁ ≫ δ₀ = 𝟙 Y' := by
    refine (hα' Y').1 (Subtype.ext (congrArg₂ Prod.mk ?_ ?_))
    · show (δ₁ ≫ δ₀) ≫ α' = 𝟙 Y' ≫ α'
      calc (δ₁ ≫ δ₀) ≫ α' = δ₁ ≫ δ₀ ≫ α' := Category.assoc _ _ _
        _ = δ₁ ≫ α := congrArg (fun t => δ₁ ≫ t) hδ₀1
        _ = α' := hδ₁1
        _ = 𝟙 Y' ≫ α' := (Category.id_comp _).symm
    · show (unTrPre P Fc).Base (δ₁ ≫ δ₀) = (unTrPre P Fc).Base (𝟙 Y')
      calc (unTrPre P Fc).Base (δ₁ ≫ δ₀)
          = (unTrPre P Fc).Base δ₁ ≫ (unTrPre P Fc).Base δ₀ := (unTrPre P Fc).Base_comp _ _
        _ = @inv _ _ _ _ hh hhiso ≫ hh := congrArg₂ (fun x y => x ≫ y) hδ₁2 hδ₀2
        _ = 𝟙 _ := IsIso.inv_hom_id hh
        _ = (unTrPre P Fc).Base (𝟙 Y') := ((unTrPre P Fc).Base_id _).symm
  -- ★手 5: `ε₀ ≫ β' = β ≫ δ₀`
  have hmid : ε₀ ≫ β' = β ≫ δ₀ := by
    refine (hα' X).1 (Subtype.ext (congrArg₂ Prod.mk ?_ ?_))
    · show (ε₀ ≫ β') ≫ α' = (β ≫ δ₀) ≫ α'
      calc (ε₀ ≫ β') ≫ α' = β ≫ α := hcancel.symm
        _ = β ≫ δ₀ ≫ α' := congrArg (fun t => β ≫ t) hδ₀1.symm
        _ = (β ≫ δ₀) ≫ α' := (Category.assoc _ _ _).symm
    · show (unTrPre P Fc).Base (ε₀ ≫ β') = (unTrPre P Fc).Base (β ≫ δ₀)
      calc (unTrPre P Fc).Base (ε₀ ≫ β')
          = (unTrPre P Fc).Base β ≫ hh := by
            rw [hhdef, ← Category.assoc, IsIso.hom_inv_id, Category.id_comp]
        _ = (unTrPre P Fc).Base β ≫ (unTrPre P Fc).Base δ₀ :=
            congrArg (fun t => (unTrPre P Fc).Base β ≫ t) hδ₀2.symm
        _ = (unTrPre P Fc).Base (β ≫ δ₀) := ((unTrPre P Fc).Base_comp _ _).symm
  -- ★組み立て
  refine ⟨⟨δ₀, δ₁, hδid, hδid'⟩, asIso ε₀, hδ₁1.symm, ?_, hε₀.symm⟩
  show β' = inv ε₀ ≫ β ≫ δ₀
  rw [← hmid, ← Category.assoc, IsIso.inv_hom_id, Category.id_comp]

include P in
/-- ★★**(i)(c)** `(𝒞^pl-bk)_A → 𝒟_{A_𝒟}` は圏同値。

★★**3 つとも pull-back 射の全単射性から直接出る**:
- **忠実性**: `Over.w` で `f ≫ W.hom = Z.hom` が言えるので、
  `Base` が一致すれば `IsPullBack W.hom` の**単射性**で消える
- **充満性**: `IsPullBack W.hom` の**全射性**で `f` を作り、
  ★**在庫の `isPullBack_of_comp_left`**(pull-back の左簡約)で `f` 自身が pull-back と分かる
- **本質的全射性**: `istr_plBkEquiv` の証人を ★`unTr_isPullBack_of_istr` で運ぶ -/
theorem unTr_plBkEquiv (Fc : FrobenioidCore P) :
    ∀ A : UnTr P, (plBkOverFunctor (unTrPre P Fc) A).IsEquivalence := by
  intro A
  haveI hfaith : (plBkOverFunctor (unTrPre P Fc) A).Faithful := by
    constructor
    intro Z W f g hfg
    have hb : (unTrPre P Fc).Base f.left.hom = (unTrPre P Fc).Base g.left.hom :=
      congrArg CommaMorphism.left hfg
    have hwf : f.left.hom ≫ W.hom.hom = Z.hom.hom :=
      congrArg InducedWideCategory.Hom.hom (Over.w f)
    have hwg : g.left.hom ≫ W.hom.hom = Z.hom.hom :=
      congrArg InducedWideCategory.Hom.hom (Over.w g)
    have h2 : f.left.hom = g.left.hom :=
      (W.hom.property Z.left.obj).1
        (Subtype.ext (congrArg₂ Prod.mk (hwf.trans hwg.symm) hb))
    exact Over.OverMorphism.ext (InducedWideCategory.Hom.ext h2)
  haveI hfull : (plBkOverFunctor (unTrPre P Fc) A).Full := by
    constructor
    intro Z W h
    have hcond : (unTrPre P Fc).Base Z.hom.hom
        = h.left ≫ (unTrPre P Fc).Base W.hom.hom := (Over.w h).symm
    obtain ⟨f, hf⟩ := (W.hom.property Z.left.obj).2 ⟨(Z.hom.hom, h.left), hcond⟩
    have hfp : (f ≫ W.hom.hom, (unTrPre P Fc).Base f) = (Z.hom.hom, h.left) :=
      Subtype.ext_iff.mp hf
    have hf1 : f ≫ W.hom.hom = Z.hom.hom := congrArg Prod.fst hfp
    have hf2 : (unTrPre P Fc).Base f = h.left := congrArg Prod.snd hfp
    have hfpb : IsPullBack (unTrPre P Fc) f :=
      isPullBack_of_comp_left (unTrPre P Fc) f W.hom.hom W.hom.property
        (by rw [hf1]; exact Z.hom.property)
    exact ⟨Over.homMk (⟨f, hfpb⟩ : Z.left ⟶ W.left) (InducedWideCategory.Hom.ext hf1),
      Over.OverMorphism.ext hf2⟩
  haveI hess : (plBkOverFunctor (unTrPre P Fc) A).EssSurj := by
    constructor
    intro Y
    haveI := istr_plBkEquiv P Fc (show Istr P from A)
    obtain ⟨Z, hZ⟩ : ∃ Z : Over (⟨show Istr P from A⟩ : PlBk (istrPre P Fc)),
        Z = (plBkOverFunctor (istrPre P Fc) (show Istr P from A)).objPreimage Y := ⟨_, rfl⟩
    have hiZ : (plBkOverFunctor (istrPre P Fc) (show Istr P from A)).obj Z ≅ Y := by
      rw [hZ]
      exact (plBkOverFunctor (istrPre P Fc) (show Istr P from A)).objObjPreimageIso Y
    refine ⟨Over.mk (⟨(istrToUnTr P).map Z.hom.hom,
      unTr_isPullBack_of_istr P Fc Z.hom.hom Z.hom.property⟩ :
        (⟨show UnTr P from Z.left.obj⟩ : PlBk (unTrPre P Fc))
          ⟶ (⟨A⟩ : PlBk (unTrPre P Fc))), ?_⟩
    exact ⟨hiZ⟩
  exact ⟨hfaith, hfull, hess⟩

/-! ## ★8. ★★★★**`𝒞^un-tr` は Frobenioid の core 21 条をすべて満たす**

原文 (FrdI p.60):
> which is faithful and essentially surjective; moreover, this functor determines
-/

include P in
/-- ★★★★**[FrdI] Proposition 3.3, (iv)** —— `𝒞^un-tr` の `FrobenioidCore` 21 条。 -/
theorem unTr_frobenioidCore (Fc : FrobenioidCore P) : FrobenioidCore (unTrPre P Fc) where
  baseSurj := unTr_baseSurj P Fc
  preStepSpan := unTr_preStepSpan P Fc
  plBkEquiv := unTr_plBkEquiv P Fc
  frobDegSurj := unTr_frobDegSurj P Fc
  frobDegUniq := unTr_frobDegUniq P Fc
  coAngularComp := unTr_coAngularComp P Fc
  coAngularOfPreStep := unTr_coAngularOfPreStep P Fc
  otriFwd := unTr_otriFwd P Fc
  otriBwd := unTr_otriBwd P Fc
  otriBase := unTr_otriBase P Fc
  arbFactor := unTr_arbFactor P Fc
  arbFactorUniq := unTr_arbFactorUniq P Fc
  pullBackLB := unTr_pullBackLB P Fc
  preStepMono := unTr_preStepMono P Fc
  preStepFactor := unTr_preStepFactor P Fc
  preStepFactorUniq := unTr_preStepFactorUniq P Fc
  preStepFactor' := unTr_preStepFactor' P Fc
  preStepFactorUniq' := unTr_preStepFactorUniq' P Fc
  faithfulUpToUnits := unTr_faithfulUpToUnits P Fc
  isotropicHullExists := unTr_isotropicHullExists P Fc
  isotropicClosed := unTr_isotropicClosed P Fc

/-! ## ★8-b. `Frobenioid` の残り 2 フィールド

★★`𝒞^un-tr` では co-angular が自動なので、`𝒞^coa-pre` は
**pre-step の広部分圏**に潰れる。★行き先 `Order(Φ(A))` は `𝒞`・`𝒞^istr`・`𝒞^un-tr`
で**同じ圏**である(底も `Φ` も変わらないから)。 -/

section CoaPre

variable (Fc : FrobenioidCore P)

/-- ★co-angular pre-step の広部分圏が乗法的であること(`letI` 用)。 -/
theorem unTr_coaPreProp_isMultiplicative :
    MorphismProperty.IsMultiplicative (coaPreProp (unTrPre P Fc)) :=
  coaPreProp_isMultiplicative (unTrPre P Fc) (unTr_frobenioidCore P Fc).coAngularComp

end CoaPre

/-! ### ★残る 2 フィールドの設計(2026-08-17、紙の上で確認)

★★**行き先 `OrderCat (Φ.val (A_𝒟))` は `𝒞`・`𝒞^istr`・`𝒞^un-tr` で同じ圏**である
(底も `Φ` も商で変わらない)。★したがって 3 つの性質は次のように取れる:

| 性質 | 手 |
|---|---|
| **忠実性** | ★**直接** —— `Under.w` の 2 式を、`Z.hom` が epi(`unTr_totEpi`)なので消す |
| **充満性** | ★`Z.hom.hom` / `W.hom.hom` を `𝒞^istr` へ持ち上げ(`istrToUnTr` は full)、`istr_coaPreUnderEquiv` の充満性で射を作り、押し戻す。★**`Div` は商で不変**なので `Order` 側の対象が一致する |
| **本質的全射性** | ★`istr_coaPreUnderEquiv` の証人を **`istrToUnTr` で押し出す**(co-angular pre-step は商を渡る) |

★★**`coaPreOverEquiv` も同型**だが、`inv (Base ...)` が入るぶん配管が重い。

★★**規律どおり、通っていないものは置いていない。** -/

/-! ## ★9. 実装の記録(2026-08-17)

★`FrobenioidCore (unTrPre)` の 21 フィールドのうち **18 が埋まった**。
★残る 3 つ **`plBkEquiv` / `arbFactorUniq` / `pullBackLB`** は、
★★**どれも同じ 1 本**に帰着する:

```
★★「`𝒞^un-tr` の pull-back 射は `𝒞^istr` の pull-back 射から来る」
   IsPullBack (unTrPre P Fc) f → ∃ α₀, (istrToUnTr P).map α₀ = f ∧ IsPullBack (istrPre P Fc) α₀
```

★これは原文 (iv) の最後の文 **「if and only if it **arises from** such an arrow of
`𝒞^istr`」の pull-back についての `⟹` の側**そのものである。
★逆(`⟸`)は上の `unTr_isPullBack_of_istr` で取れている。

### ★なぜ難しいか(測定)

★`𝒞^istr → 𝒞` の同じ向きの移送(`istr_isPullBack_to`、`Prop19.lean`)は
★**isotropification が包含関手の左随伴であること**(`hullHomEquiv`)で埋めていた
——「すべての `Z`」と「isotropic な `Z`」の差を随伴が吸収する。

★★**しかし `𝒞^istr → 𝒞^un-tr` の差は「商」であって随伴ではない。**
`Hom_{un-tr}(X,A) = Hom_{istr}(X,A) / ≈` なので、
★`𝒞^un-tr` の全単射性から `𝒞^istr` の全単射性は**そのままでは出ない**
(単射性を戻すと「等しい」ではなく「unit-equivalent」しか言えない)。

### ★試して**駄目だった**筋(記録)

1. ★**「どの持ち上げでも pull-back」に強める** —— ✗。
   `IsPullBack (istrPre)` は unit-equivalence で**不変ではない**:
   `f₁ ≫ α₀ ≈ f₂ ≫ α₀` から `f₁ ≫ α₀ = f₂ ≫ α₀` は出ない。
2. ★**`𝒞^un-tr` の中で直接 `pullBackLB` を示す** —— ✗(未達)。
   `Definition 1.3, (iv)(b)` は原文では**公理**であり、
   全単射性から `degFr = 1` を出すには「次数 1 の射が像に居る」ことが要るが、
   `X` の取り方でそれを作る手が見つからなかった。

### ★次の一手の候補(未検証)

- ★`arbFactor` で `α` を分解し、**pull-back 部分だけを持ち上げる**
- ★`Proposition 3.3, (ii)`(unit-equivalent ⟺ `𝔽_Φ` で同じ)を使って、
  **持ち上げを単元で補正する**(`f₀ ≫ α₀ ≈ g₀` を等式にする)
-/

end ABC3.Found.FrdI
