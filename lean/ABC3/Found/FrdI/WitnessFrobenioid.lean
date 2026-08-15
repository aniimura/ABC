import ABC3.Found.FrdI.Frobenioid
import ABC3.Found.FrdI.Witness

/-!
# ★`wP` は Frobenioid か —— [FrdI] Definition 1.3 の非退化判定

`Witness.lean` で組み立てた本物の pre-Frobenioid `wP`
(`𝒟 = Vee`、`Φ = ℕ` への定数関手、`𝒞 = 𝔽_Φ`)が
`Definition 1.3` の条件を満たすかを、**条件ごとに**判定する。

## ★このモデルの構造(先に測る)

`wC` の射 `A ⟶ B` は3つ組 `(base, div, deg) ∈ Hom_Vee(A.base, B.base) × ℕ × ℕ≥1` で、
`Vee` の hom は subsingleton なので、**実質は `(div, deg) ∈ ℕ × ℕ≥1`** である。
合成は `(b,d,n) ≫ (b',d',n') = (b ≫ b', d' + n'·d, n'·n)`。

この構造から3つの事実が出る(下で証明する):

1. ★**isometric な pre-step はすべて同型**(`Witness.lean` の
   `wIsIso_of_isometric_preStep`)
2. ★したがって **`wC` のすべての射が co-angular** —— co-angular の条件は
   「`β` が同型」を要求するが、`β` は isometric pre-step なので自動的に同型
3. ★**`IsPullBack φ ⟺ φ.div = 0 ∧ φ.deg = 1`**

## ★実装上の注意(2度目に踏んだ壁)

`wΦ.val A` は `ℕ` と**定義的に等しい**が、インスタンス解決は `reducible` 透明度しか
使わないので `HAdd (wΦ.val A) ℕ` が見つからない。そこで零因子・次数を
**`ℕ` / `ℕ+` として読み出す関数** `wd` / `wn` を通す。
これらは `rfl` で定義でき、以後の等式はすべて `ℕ` の中の話になる。
-/

namespace ABC3.Found.FrdI

open CategoryTheory Opposite

/-! ### ★`ℕ` として読み出す -/

/-- `wC` の射の零因子を `ℕ` として読む。 -/
def wd {A B : wC} (φ : A ⟶ B) : ℕ := ElemFrobCat.Hom.div φ

/-- `wC` の射の Frobenius 次数。 -/
def wn {A B : wC} (φ : A ⟶ B) : ℕ+ := ElemFrobCat.Hom.deg φ

@[simp] theorem wd_comp {A B E : wC} (ψ : A ⟶ B) (φ : B ⟶ E) :
    wd (ψ ≫ φ) = wd φ + ((wn φ : ℕ+) : ℕ) • wd ψ := rfl

@[simp] theorem wn_comp {A B E : wC} (ψ : A ⟶ B) (φ : B ⟶ E) :
    wn (ψ ≫ φ) = wn φ * wn ψ := rfl

@[simp] theorem wd_id (A : wC) : wd (𝟙 A) = 0 := rfl
@[simp] theorem wn_id (A : wC) : wn (𝟙 A) = 1 := rfl

@[simp] theorem wd_mk {A B : wC} (b : A.base ⟶ B.base) (d : ℕ) (n : ℕ+) :
    wd (⟨b, d, n⟩ : A ⟶ B) = d := rfl

@[simp] theorem wn_mk {A B : wC} (b : A.base ⟶ B.base) (d : ℕ) (n : ℕ+) :
    wn (⟨b, d, n⟩ : A ⟶ B) = n := rfl

/-- `Vee` の hom は subsingleton なので、`wC` の射は `wd` と `wn` で決まる。 -/
theorem wHom_ext {A B : wC} {f g : A ⟶ B} (hd : wd f = wd g) (hn : wn f = wn g) : f = g :=
  ElemFrobCat.Hom.ext (Subsingleton.elim _ _) hd hn

theorem wIsLinear_iff {A B : wC} (φ : A ⟶ B) : IsLinear wP φ ↔ wn φ = 1 := Iff.rfl
theorem wIsIsometric_iff {A B : wC} (φ : A ⟶ B) : IsIsometric wP φ ↔ wd φ = 0 := Iff.rfl

/-! ### ★基本構造 -/

/-- ★**`wC` のすべての射は co-angular**。

co-angular の条件は「分解 `φ = γ ≫ β ≫ α` で `β` が isometric pre-step なら
`β` は同型」だが、`wIsIso_of_isometric_preStep` によりそれは**無条件に成り立つ**。 -/
theorem wIsCoAngular {A B : wC} (φ : A ⟶ B) : IsCoAngular wP φ :=
  fun _ _ _ β _ _ _ hiso hstep _ => wIsIso_of_isometric_preStep β hiso hstep

/-- ★**`wd = 0` かつ `wn = 1` なら pull-back morphism**。 -/
theorem wIsPullBack_of {A B : wC} (φ : A ⟶ B) (hd : wd φ = 0) (hn : wn φ = 1) :
    IsPullBack wP φ := by
  intro X
  constructor
  · intro f₁ f₂ h
    have h1 : f₁ ≫ φ = f₂ ≫ φ :=
      congrArg (fun p => (p : (X ⟶ B) × _).1) (congrArg Subtype.val h)
    refine wHom_ext ?_ ?_
    · have h2 : wd (f₁ ≫ φ) = wd (f₂ ≫ φ) := congrArg wd h1
      simpa [hd, hn] using h2
    · have h2 : wn (f₁ ≫ φ) = wn (f₂ ≫ φ) := congrArg wn h1
      simpa [hn] using h2
  · rintro ⟨⟨g, h'⟩, -⟩
    refine ⟨⟨h', wd g, wn g⟩, Subtype.ext (Prod.ext ?_ rfl)⟩
    refine wHom_ext ?_ ?_
    · show wd (_ ≫ φ) = wd g
      simp [hd, hn]
    · show wn (_ ≫ φ) = wn g
      simp [hn]

/-- ★**逆向き** —— pull-back morphism なら `wd = 0` かつ `wn = 1`。

`X := A` で全射性を使う。`(⟨φ.base, 0, 1⟩, 𝟙)` の原像 `f` を取ると
`wd φ + wn φ · wd f = 0` と `wn φ * wn f = 1` が出る。 -/
theorem wPullBack_wd_wn {A B : wC} (φ : A ⟶ B) (h : IsPullBack wP φ) :
    wd φ = 0 ∧ wn φ = 1 := by
  obtain ⟨f, hf⟩ := (h A).2
    ⟨(⟨ElemFrobCat.Hom.base φ, (0 : ℕ), 1⟩, 𝟙 _), Subsingleton.elim _ _⟩
  have h1 : f ≫ φ = (⟨ElemFrobCat.Hom.base φ, (0 : ℕ), 1⟩ : A ⟶ B) :=
    congrArg (fun p => (p : (A ⟶ B) × _).1) (congrArg Subtype.val hf)
  have hdiv : wd φ + ((wn φ : ℕ+) : ℕ) * wd f = 0 := by
    have h2 : wd (f ≫ φ) = wd (⟨ElemFrobCat.Hom.base φ, (0 : ℕ), 1⟩ : A ⟶ B) := congrArg wd h1
    simpa [smul_eq_mul] using h2
  have hdeg : wn φ * wn f = 1 := by
    have h2 : wn (f ≫ φ) = wn (⟨ElemFrobCat.Hom.base φ, (0 : ℕ), 1⟩ : A ⟶ B) := congrArg wn h1
    simpa using h2
  refine ⟨by omega, ?_⟩
  have hc := congrArg (fun n : ℕ+ => (n : ℕ)) hdeg
  simp only [PNat.mul_coe, PNat.one_coe] at hc
  exact PNat.coe_eq_one_iff.mp (Nat.eq_one_of_mul_eq_one_right hc)

/-- ★**`IsPullBack` の完全な特徴づけ**(このモデルでの実測)。

★これは「`wP` は本当に何かを区別している」ことの証拠でもある——
pull-back morphism は全射ではなく、`wd`・`wn` で切り出される真部分クラスである。 -/
theorem wIsPullBack_iff {A B : wC} (φ : A ⟶ B) :
    IsPullBack wP φ ↔ wd φ = 0 ∧ wn φ = 1 :=
  ⟨wPullBack_wd_wn φ, fun h => wIsPullBack_of φ h.1 h.2⟩

/-- ★`wDeg2`(次数 2)は **pull-back morphism でない** —— 上の特徴づけの非退化。 -/
theorem wNot_pullBack_deg2 : ¬ IsPullBack wP wDeg2 := by
  intro h
  have h2 : (2 : ℕ+) = 1 := (wPullBack_wd_wn _ h).2
  exact absurd h2 (by decide)

/-! ### ★`Definition 1.3` の条件を1つずつ判定する -/

/-- **(i)(a)** `Vee` の各対象は Frobenius-trivial な対象の底として現れる。

`ζ n := (𝟙, 0, n)` がモノイド準同型 `ℕ≥1 → End A` を与える。 -/
theorem wBaseSurj (Y : Vee) :
    ∃ A : wC, IsFrobeniusTrivial wP A ∧ Nonempty ((wP.toElem.obj A).base ≅ Y) := by
  refine ⟨⟨Y⟩, ⟨⟨⟨fun n => ⟨𝟙 Y, (0 : ℕ), n⟩, rfl⟩, fun m n => ?_⟩, fun n => rfl,
    fun n => ⟨rfl, ⟨wIsCoAngular _, rfl⟩, ?_⟩⟩, ⟨Iso.refl _⟩⟩
  · refine wHom_ext ?_ ?_
    · show (0 : ℕ)
        = wd ((⟨𝟙 Y, (0 : ℕ), n⟩ : (⟨Y⟩ : wC) ⟶ ⟨Y⟩) ≫ ⟨𝟙 Y, (0 : ℕ), m⟩)
      simp
    · show m * n
        = wn ((⟨𝟙 Y, (0 : ℕ), n⟩ : (⟨Y⟩ : wC) ⟶ ⟨Y⟩) ≫ ⟨𝟙 Y, (0 : ℕ), m⟩)
      simp
  · show IsIso (𝟙 Y); infer_instance

/-- **(ii)** 各次数の Frobenius 型射が存在する。 -/
theorem wFrobDegSurj (A : wC) (n : ℕ+) :
    ∃ (B : wC) (φ : A ⟶ B), IsFrobeniusType wP φ ∧ wP.degFr φ = n := by
  refine ⟨A, ⟨𝟙 _, (0 : ℕ), n⟩, ⟨⟨wIsCoAngular _, rfl⟩, ?_⟩, rfl⟩
  show IsIso (𝟙 A.base)
  infer_instance

/-- **(iii)(a)** co-angular 射は合成で閉じる(このモデルでは全射が co-angular)。 -/
theorem wCoAngularComp {A B E : wC} (ψ : A ⟶ B) (φ : B ⟶ E) :
    IsCoAngular wP ψ → IsCoAngular wP φ → IsCoAngular wP (ψ ≫ φ) :=
  fun _ _ => wIsCoAngular _

/-- **(iii)(b)** -/
theorem wCoAngularOfPreStep {A' A : wC} (α : A' ⟶ A) :
    IsCoAngular wP α → IsPreStep wP α → ∀ φ : A' ⟶ A, IsCoAngular wP φ :=
  fun _ _ φ => wIsCoAngular φ

/-- **(iv)(a)** 任意射の分解 `φ = γ ≫ β ≫ α`。

`φ = (b, d, n)` に対し `γ := (𝟙, 0, n)`(Frobenius 型)、`β := (𝟙, d, 1)`(pre-step)、
`α := (b, 0, 1)`(pull-back)と取る。 -/
theorem wArbFactor {A B : wC} (φ : A ⟶ B) :
    ∃ (X Y : wC) (γ : A ⟶ X) (β : X ⟶ Y) (α : Y ⟶ B),
      φ = γ ≫ β ≫ α ∧ IsFrobeniusType wP γ ∧ IsPreStep wP β ∧ IsPullBack wP α := by
  refine ⟨A, A, ⟨𝟙 _, (0 : ℕ), wn φ⟩, ⟨𝟙 _, wd φ, 1⟩, ⟨ElemFrobCat.Hom.base φ, (0 : ℕ), 1⟩,
    wHom_ext (by simp) (by simp), ⟨⟨wIsCoAngular _, rfl⟩, ?_⟩, ⟨rfl, ?_⟩,
    wIsPullBack_of _ rfl rfl⟩
  · show IsIso (𝟙 A.base); infer_instance
  · show IsIso (𝟙 A.base); infer_instance

/-- **(iv)(b)** pull-back morphism は LB-invertible かつ linear。 -/
theorem wPullBackLB {A B : wC} (α : A ⟶ B) (h : IsPullBack wP α) :
    IsLBInvertible wP α ∧ IsLinear wP α :=
  ⟨⟨wIsCoAngular α, (wPullBack_wd_wn α h).1⟩, (wPullBack_wd_wn α h).2⟩

/-- **(v)(a)** pre-step は monomorphism(このモデルでは**すべての射**が mono)。 -/
theorem wMono {A B : wC} (φ : A ⟶ B) : Mono φ := by
  refine ⟨fun {Z} f g h => ?_⟩
  have hd : wd (f ≫ φ) = wd (g ≫ φ) := congrArg wd h
  have hn : wn (f ≫ φ) = wn (g ≫ φ) := congrArg wn h
  simp only [wd_comp, smul_eq_mul] at hd
  simp only [wn_comp] at hn
  refine wHom_ext (Nat.eq_of_mul_eq_mul_left (wn φ).pos (Nat.add_left_cancel hd)) ?_
  exact mul_left_cancel hn

/-- **(v)(b)** pre-step は「co-angular pre-step ∘ isometric pre-step」に分解する。 -/
theorem wPreStepFactor {A B : wC} (φ : A ⟶ B) (hφ : IsPreStep wP φ) :
    ∃ (X : wC) (β : A ⟶ X) (α : X ⟶ B),
      φ = β ≫ α ∧ IsCoAngular wP β ∧ IsPreStep wP β ∧ IsIsometric wP α ∧ IsPreStep wP α := by
  refine ⟨A, ⟨𝟙 _, wd φ, 1⟩, ⟨ElemFrobCat.Hom.base φ, (0 : ℕ), 1⟩,
    wHom_ext (by simp) (by simp [show wn φ = 1 from hφ.1]), wIsCoAngular _, ⟨rfl, ?_⟩,
    rfl, ⟨rfl, hφ.2⟩⟩
  show IsIso (𝟙 A.base); infer_instance

/-- **(v)(c)** pre-step は「isometric pre-step ∘ co-angular pre-step」にも分解する。 -/
theorem wPreStepFactor' {A B : wC} (φ : A ⟶ B) (hφ : IsPreStep wP φ) :
    ∃ (X : wC) (β : A ⟶ X) (α : X ⟶ B),
      φ = β ≫ α ∧ IsIsometric wP β ∧ IsPreStep wP β ∧ IsCoAngular wP α ∧ IsPreStep wP α := by
  refine ⟨B, ⟨ElemFrobCat.Hom.base φ, (0 : ℕ), 1⟩, ⟨𝟙 _, wd φ, 1⟩,
    wHom_ext (by simp) (by simp [show wn φ = 1 from hφ.1]), rfl, ⟨rfl, hφ.2⟩,
    wIsCoAngular _, ⟨rfl, ?_⟩⟩
  show IsIso (𝟙 B.base); infer_instance

/-- **(vii)(b)** isotropic な対象からの射の終域も isotropic(このモデルでは全対象)。 -/
theorem wIsotropicClosed {A B : wC} (φ : A ⟶ B) : IsIsotropic wP A → IsIsotropic wP B :=
  fun _ => wIsotropic B

/-- **(vii)(a)** isotropic hull が存在する —— `𝟙 A` がそれである。 -/
theorem wIsotropicHullExists (A : wC) : ∃ (B : wC) (φ : A ⟶ B), IsIsotropicHull wP φ := by
  refine ⟨A, 𝟙 A, isIsometric_id wP A, isPreStep_id wP A, wIsotropic A, ?_⟩
  intro Cc _ γ
  exact ⟨γ, (Category.id_comp γ).symm, fun y hy => by simpa using hy.symm⟩

end ABC3.Found.FrdI
