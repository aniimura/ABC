import ABC3.Found.FrdI.Remark311
import ABC3.Found.FrdI.Witness
import Mathlib.CategoryTheory.SingleObj

/-!
# [FrdI] Example 3.10 —— Non-slim Base Categories

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.72。

原文 (FrdI p.72):
> Non-slim Base Categories. Let G be a group, whose center

## ★この例の中身(測定)

**主張は 6**:

| # | 内容 |
|---|---|
| 1 | `𝒞` の自己射モノイドは `G × 𝔽` |
| 2 | `Φ`(定数 `ℤ≥0`)は non-dilating |
| 3 | `𝒞` は Frobenius-normalized かつ isotropic 型で、**group-like 型ではない** |
| 4 | `𝒟` は FSM-type、したがって FSMFF-type |
| 5 | ★**`𝒞` は standard 型** |
| 6 | ★★**`α : 𝔽 → Z(G)` が非自明なら、`(g,f) ↦ (g·α(f), f)` が定める自己同値は
base-identity な Frobenius 型自己射を保たない** |

★**6 が眼目**である —— `Theorem 3.4, (v)` が `𝒟` の Frobenius-slim 性を
要求することの必然性を示す反例になっている。

## ★合成の向き

`SingleObj G` では `f ≫ g = g * f`(mathlib の規約)。
★`End A` の積は `x * y = y ≫ x` なので、`base` 成分では `G` の積とそのまま合う。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

variable (G : Type) [Group G]

/-! ## ★底圏 `𝒟 = SingleObj G` -/

/-- ★**底圏** —— 唯一の対象の自己射モノイドが `G` である 1 対象圏。 -/
abbrev D310 : Type := SingleObj G

instance : Groupoid (D310 G) := inferInstanceAs (Groupoid (SingleObj G))

instance : IsConnected (D310 G) := by
  refine IsConnected.of_induct (j₀ := (SingleObj.star G)) ?_
  intro p hp₀ _ A
  obtain rfl : A = SingleObj.star G := Subsingleton.elim _ _
  exact hp₀

/-- ★**`𝒟` の射はすべて同型** —— `G` は群だから。 -/
theorem d310_isIso {A B : D310 G} (f : A ⟶ B) : IsIso f := inferInstance

/-- ★**`𝒟` は totally epimorphic** —— 同型は epi。 -/
theorem d310_totEpi : IsTotallyEpimorphic (D310 G) := by
  intro A B f
  haveI := d310_isIso G f
  infer_instance

/-- ★★**主張 4** —— `𝒟` は **FSM-type**(射がすべて同型だから)。 -/
theorem d310_fsmType : IsOfFSMType (D310 G) := fun _ _ β _ => d310_isIso G β

/-- ★★**主張 4** —— したがって **FSMFF-type**。 -/
theorem d310_fsmffType : IsOfFSMFFType (D310 G) :=
  isOfFSMFFType_of_isOfFSMType (d310_fsmType G)

/-! ## ★`Φ` と `𝒞 = 𝔽_Φ` -/

/-- ★**`Φ`** —— `𝒟` の唯一の対象に `ℤ≥0 = ℕ` を、どの射にも恒等自己同型を割り当てる。 -/
abbrev Phi310 : MonoidOn.{0, 0, 0} (D310 G) := MonoidOn.const (D310 G) ℕ

/-- ★★**主張 2** —— `Φ` は **non-dilating**(定数関手なので誘導射は恒等)。 -/
theorem phi310_nonDilating : MonoidOn.IsNonDilatingOn (Phi310 G) := by
  intro A e
  have h : (Phi310 G).map e = AddMonoidHom.id ℕ := by ext x; rfl
  rw [h]
  exact isNonDilating_id (M := ℕ)

/-- ★★**主張 2**(`Φ^char` 側)—— `P310` が担ぐモノイドは `Φ^char` なのでこちらが要る。 -/
theorem phi310_charOn_nonDilating : MonoidOn.IsNonDilatingOn (Phi310 G).charOn := by
  intro A e
  have h : (Phi310 G).charOn.map e = AddMonoidHom.id (MChar ℕ) := by
    ext x
    obtain ⟨m, rfl⟩ := toChar_surjective ℕ x
    rfl
  rw [h]
  exact isNonDilating_id (M := MChar ℕ)

/-- ★**`𝒞 = 𝔽_Φ` の pre-Frobenioid 構造**。 -/
abbrev P310 : PreFrobenioid (ElemFrobCat (Phi310 G)) (Phi310 G).charOn :=
  elemPreFrobenioid (Phi310 G) (d310_totEpi G) (fun _ => isPreDivisorial_nat)

/-! ## ★主張 1 —— 自己射モノイドは `G × 𝔽` -/

/-- ★★**主張 1** —— `𝒞` の自己射モノイドは `G × 𝔽`。

★`End` の積は `x * y = y ≫ x`、`SingleObj G` では `f ≫ g = g * f` なので、
`base` 成分は `G` の積とそのまま合う。`div`・`deg` 成分は `𝔽 = ElemFrob ℕ` の積。 -/
def ex310EndEquiv (A : ElemFrobCat (Phi310 G)) :
    End A ≃* G × ElemFrob ℕ where
  toFun f := (f.base, ⟨f.div, f.deg⟩)
  invFun p := ⟨p.1, p.2.div, p.2.deg⟩
  left_inv _ := rfl
  right_inv _ := rfl
  map_mul' x y := rfl

/-! ## ★主張 3 —— 型の判定 -/

/-- ★★**主張 3** —— `𝒞` は **isotropic 型**。 -/
theorem ex310_isotropicType : IsOfIsotropicType (P310 G) :=
  fun A => elemFrob_isotropic (Phi310 G) (d310_totEpi G) (fun _ => isPreDivisorial_nat) A

/-- ★★**主張 3** —— `𝒞` は **Frobenius-normalized 型**。 -/
theorem ex310_frobNormalizedType : IsOfFrobeniusNormalizedType (P310 G) :=
  fun A => elemFrob_frobeniusNormalized (Phi310 G) (d310_totEpi G)
    (fun _ => isPreDivisorial_nat) A

/-- ★`ℕ^char` も group-like でない —— `M^char` は sharp なので、
すべてが可逆なら全部 `0` になり `ℕ` が group-like になってしまう。 -/
theorem not_isGroupLike_mChar_nat : ¬ IsGroupLike (MChar ℕ) := by
  intro h
  refine not_isGroupLike_nat ?_
  intro x
  exact isSharp_mChar ℕ x ((isGroupLike_iff (MChar ℕ)).mp h x)

/-- ★★**主張 3** —— `𝒞` は **group-like 型ではない**。 -/
theorem ex310_not_groupLikeType : ¬ IsOfGroupLikeType (P310 G) := by
  intro h
  exact not_isGroupLike_mChar_nat (h ⟨SingleObj.star G⟩)

/-! ## ★主張 5 —— `𝒞` は standard 型 -/

/-- ★★★**主張 5** —— `𝒞` は **standard 型**。

★(a) isotropic 型なので quasi-isotropic 型(`Remark 3.1.1`)、
かつ Frobenius-isotropic 型(恒等射が Frobenius 型)。
★(b) group-like 型でないので**前件が偽**、含意は自明に真。
★(c)(d)(e) は上で取った。 -/
theorem ex310_frobenioid : Frobenioid (P310 G) :=
  elemFrob_isFrobenioid (Phi310 G) (d310_totEpi G) (fun _ => isPreDivisorial_nat)

/-- ★`𝒞` の `FrobenioidCore`。 -/
theorem ex310_core : FrobenioidCore (P310 G) :=
  elemFrob_frobenioidCore (Phi310 G) (d310_totEpi G) (fun _ => isPreDivisorial_nat)

theorem ex310_standardType :
    IsOfStandardType (D310 G) (ElemFrobCat (Phi310 G)) (P310 G) (ex310_core G) where
  quasiIsotropic :=
    isOfQuasiIsotropicType_of_isOfIsotropicType (P310 G) (ex310_core G) (ex310_isotropicType G)
  frobIsotropic := fun A =>
    ⟨A, 𝟙 A, isFrobeniusType_of_isIso (P310 G) (𝟙 A), ex310_isotropicType G A⟩
  groupLikeCompact := fun hgl => absurd hgl (ex310_not_groupLikeType G)
  frobNormalized := ex310_frobNormalizedType G
  baseFSMFF := d310_fsmffType G
  phiNonDilating := phi310_charOn_nonDilating G

/-! ## ★★主張 6 —— 捻れた自己同値は base-identity Frobenius 型を保たない

原文 (FrdI p.72):
> determines a self-equivalence C
-/

section Twist

variable (α : ℕ+ →* G) (hα : ∀ n : ℕ+, α n ∈ Subgroup.center G)

/-- ★★**捻れた自己関手** —— `(g, f) ↦ (g · α(deg f), f)`。

★原文の `α : 𝔽 → Z(G)` は `𝔽 ↠ ℕ≥1` を経由するので、
**次数だけに依る**捻れである。 -/
def ex310Twist : ElemFrobCat (Phi310 G) ⥤ ElemFrobCat (Phi310 G) where
  obj A := A
  map {_ _} f := ⟨α f.deg * f.base, f.div, f.deg⟩
  map_id A := by
    refine ElemFrobCat.Hom.ext ?_ rfl rfl
    show α 1 * (1 : G) = (1 : G)
    rw [map_one, one_mul]
  map_comp {A B E} f g := by
    refine ElemFrobCat.Hom.ext ?_ rfl rfl
    show α (g.deg * f.deg) * (g.base * f.base)
      = (α g.deg * g.base) * (α f.deg * f.base)
    rw [map_mul]
    have hc : α f.deg * g.base = g.base * α f.deg :=
      ((Subgroup.mem_center_iff.mp (hα f.deg)) g.base).symm
    simp only [mul_assoc]
    rw [← mul_assoc (α f.deg), hc, mul_assoc]

/-- ★**捻れた関手は充満忠実かつ本質的全射**、すなわち**自己同値**である。 -/
instance : (ex310Twist G α hα).Faithful where
  map_injective {A B} {f g} h := by
    have hdiv := congrArg ElemFrobCat.Hom.div h
    have hdeg := congrArg ElemFrobCat.Hom.deg h
    have hb := congrArg ElemFrobCat.Hom.base h
    have hdeg' : f.deg = g.deg := hdeg
    have hb' : α f.deg * f.base = α g.deg * g.base := hb
    refine ElemFrobCat.Hom.ext ?_ hdiv hdeg'
    rw [hdeg'] at hb'
    exact mul_left_cancel hb'

instance : (ex310Twist G α hα).Full where
  map_surjective {A B} h := by
    refine ⟨⟨(α h.deg)⁻¹ * h.base, h.div, h.deg⟩, ?_⟩
    refine ElemFrobCat.Hom.ext ?_ rfl rfl
    show α h.deg * ((α h.deg)⁻¹ * h.base) = h.base
    rw [← mul_assoc, mul_inv_cancel, one_mul]

instance : (ex310Twist G α hα).EssSurj where
  mem_essImage A := ⟨A, ⟨Iso.refl _⟩⟩

/-- ★★★**主張 6** —— `α` が非自明なら、捻れた自己関手は
**base-identity な Frobenius 型自己射を保たない**。

★`⟨𝟙, 0, n⟩` は base-identity な Frobenius 型自己射だが、
その像の底は `α n` であって `𝟙` ではない。 -/
theorem ex310Twist_not_preserves_baseIdentity
    {n : ℕ+} (hn : α n ≠ 1) (A : ElemFrobCat (Phi310 G)) :
    ∃ φ : End A, IsBaseIdentity (P310 G) φ ∧ IsFrobeniusType (P310 G) φ ∧
      ¬ IsBaseIdentity (P310 G) ((ex310Twist G α hα).map φ) := by
  refine ⟨⟨𝟙 _, 0, n⟩, rfl, ?_, ?_⟩
  · refine (elemFrob_frobeniusType_iff (Phi310 G) (d310_totEpi G)
      (fun _ => isPreDivisorial_nat) _).mpr ⟨?_, ?_⟩
    · show toChar (0 : ℕ) = 0
      exact map_zero _
    · show IsIso (𝟙 _)
      infer_instance
  · intro h
    refine hn ?_
    have h' : α n * (1 : G) = (1 : G) := h
    rwa [mul_one] at h'

end Twist

/-! ## ★★★出典の紐付け(`.src`) -/

/-- ★★★**[FrdI] Example 3.10** —— 6 主張すべてが実装された。

| # | 主張 | 実装 |
|---|---|---|
| 1 | 自己射モノイドは `G × 𝔽` | `ex310EndEquiv` |
| 2 | `Φ` は non-dilating | `phi310_nonDilating` / `phi310_charOn_nonDilating` |
| 3 | Frobenius-normalized・isotropic 型、group-like 型でない | `ex310_frobNormalizedType` / `ex310_isotropicType` / `ex310_not_groupLikeType` |
| 4 | `𝒟` は FSM-type ゆえ FSMFF-type | `d310_fsmType` / `d310_fsmffType` |
| 5 | `𝒞` は standard 型 | `ex310_standardType`(＋ `ex310_frobenioid`) |
| 6 | 捻れた自己同値は base-identity Frobenius 型を保たない | `ex310Twist`(充満忠実・本質的全射)＋ `ex310Twist_not_preserves_baseIdentity` |

★★**6 が眼目**である —— `Theorem 3.4, (v)` が `𝒟` の Frobenius-slim 性を
要求することの必然性を、この例が示している。
★`𝒟 = SingleObj G` は `Z(G)` が非自明なら slim でない。 -/
def ex310_standardType.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 72, item := "Example 3.10",
    sectionId := "frdi-example-3-10" }

end ABC3.Found.FrdI
