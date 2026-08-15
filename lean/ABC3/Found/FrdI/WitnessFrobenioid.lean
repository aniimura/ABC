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

/-! ### ★第1段 —— 一意性系と (iii)(c) -/

/-- 底の同型から `wC` の同型を作る(`wd = 0`、`wn = 1`)。

★逆射を**明示して**構成する。`asIso` だと逆射の `wd` が制御できない。 -/
noncomputable def wIso {A B : wC} (b : A.base ⟶ B.base) (hb : IsIso b) : A ≅ B where
  hom := ⟨b, (0 : ℕ), 1⟩
  inv := ⟨@inv _ _ _ _ b hb, (0 : ℕ), 1⟩
  hom_inv_id := wHom_ext (by simp) (by simp)
  inv_hom_id := wHom_ext (by simp) (by simp)

@[simp] theorem wd_wIso_hom {A B : wC} (b : A.base ⟶ B.base) (hb : IsIso b) :
    wd (wIso b hb).hom = 0 := rfl
@[simp] theorem wn_wIso_hom {A B : wC} (b : A.base ⟶ B.base) (hb : IsIso b) :
    wn (wIso b hb).hom = 1 := rfl
@[simp] theorem wd_wIso_inv {A B : wC} (b : A.base ⟶ B.base) (hb : IsIso b) :
    wd (wIso b hb).inv = 0 := rfl
@[simp] theorem wn_wIso_inv {A B : wC} (b : A.base ⟶ B.base) (hb : IsIso b) :
    wn (wIso b hb).inv = 1 := rfl

/-- **(i)(b)** 底の同型は pre-step の span で持ち上がる。

★`Vee` の hom は subsingleton なので、等式そのものは `Subsingleton.elim` で済む。 -/
theorem wPreStepSpan (A B : wC) (α : (wP.toElem.obj A).base ⟶ (wP.toElem.obj B).base)
    (hα : IsIso α) :
    ∃ (X : wC) (φ : X ⟶ A) (ψ : X ⟶ B) (hφ : IsPreStep wP φ), IsPreStep wP ψ ∧
      α = @inv _ _ _ _ (wP.Base φ) hφ.2 ≫ wP.Base ψ := by
  refine ⟨A, 𝟙 A, ⟨α, (0 : ℕ), 1⟩, isPreStep_id wP A, ⟨rfl, hα⟩, Subsingleton.elim _ _⟩

/-- **(ii)** Frobenius 型射の本質的一意性。

Frobenius 型は isometric(`wd = 0`)かつ base-isomorphism なので、
底の同型を繋いだ `β := (inv (Base φ) ≫ Base ψ, 0, 1)` が求めるもの。 -/
theorem wFrobDegUniq (A B E : wC) (φ : A ⟶ B) (ψ : A ⟶ E)
    (hφ : IsFrobeniusType wP φ) (hψ : IsFrobeniusType wP ψ) (hn : wP.degFr φ = wP.degFr ψ) :
    ∃ β : B ⟶ E, IsIso β ∧ φ ≫ β = ψ := by
  haveI hbφ : IsIso (wP.Base φ) := hφ.2
  haveI hbψ : IsIso (wP.Base ψ) := hψ.2
  refine ⟨⟨(@inv _ _ _ _ (wP.Base φ) hbφ) ≫ wP.Base ψ, (0 : ℕ), 1⟩,
    wIsIso_of_isometric_preStep _ rfl
      ⟨rfl, (inferInstance : IsIso ((@inv _ _ _ _ (wP.Base φ) hbφ) ≫ wP.Base ψ))⟩, ?_⟩
  refine wHom_ext ?_ ?_
  · have h1 : wd φ = 0 := hφ.1.2
    have h2 : wd ψ = 0 := hψ.1.2
    simp [h1, h2]
  · have : wn φ = wn ψ := hn
    simp [this]

/-- **(iii)(c)** co-angular pre-step は `𝒪^▷` の全単射を誘導する(順方向)。

`φ ≫ β` と `α ≫ φ` の `wd` を比べると `wd β = wd α` に落ちる(`wn φ = 1` だから)。 -/
theorem wOtriFwd {A B : wC} (φ : A ⟶ B) (hca : IsCoAngular wP φ) (hst : IsPreStep wP φ)
    (α : End A) (hα : α ∈ OTri wP A) :
    ∃! β : End B, β ∈ OTri wP B ∧ (φ ≫ (β : B ⟶ B)) = (α : A ⟶ A) ≫ φ := by
  have hφn : wn φ = 1 := hst.1
  have hαn : wn (α : A ⟶ A) = 1 := hα.2
  refine ⟨(⟨𝟙 _, wd (α : A ⟶ A), 1⟩ : End B), ⟨⟨rfl, rfl⟩, ?_⟩, ?_⟩
  · exact wHom_ext (by simp [hφn, hαn, add_comm]) (by simp [hφn, hαn])
  · rintro β ⟨⟨-, hβn⟩, hβ⟩
    have hd : wd (φ ≫ (β : B ⟶ B)) = wd ((α : A ⟶ A) ≫ φ) := congrArg wd hβ
    have hβn' : wn (β : B ⟶ B) = 1 := hβn
    simp only [wd_comp, hφn, hβn', PNat.one_coe, one_smul] at hd
    exact wHom_ext (by simp [hβn']; omega) (by simp [hβn'])

/-- **(iii)(c)** 逆方向。 -/
theorem wOtriBwd {A B : wC} (φ : A ⟶ B) (hca : IsCoAngular wP φ) (hst : IsPreStep wP φ)
    (β : End B) (hβ : β ∈ OTri wP B) :
    ∃! α : End A, α ∈ OTri wP A ∧ (φ ≫ (β : B ⟶ B)) = (α : A ⟶ A) ≫ φ := by
  have hφn : wn φ = 1 := hst.1
  have hβn : wn (β : B ⟶ B) = 1 := hβ.2
  refine ⟨(⟨𝟙 _, wd (β : B ⟶ B), 1⟩ : End A), ⟨⟨rfl, rfl⟩, ?_⟩, ?_⟩
  · exact wHom_ext (by simp [hφn, hβn, add_comm]) (by simp [hφn, hβn])
  · rintro α ⟨⟨-, hαn⟩, hα⟩
    have hd : wd (φ ≫ (β : B ⟶ B)) = wd ((α : A ⟶ A) ≫ φ) := congrArg wd hα
    have hαn' : wn (α : A ⟶ A) = 1 := hαn
    simp only [wd_comp, hφn, hβn, hαn', PNat.one_coe, one_smul] at hd
    exact wHom_ext (by simp [hαn']; omega) (by simp [hαn'])

/-- **(iii)(c)** 全単射は `Base(φ)` にしか依らない。

★このモデルでは `Base` は subsingleton なので条件は空だが、
**結論の側は非自明**(`wd` の等式が `φ` に依らないことを示す必要がある)。 -/
theorem wOtriBase {A B : wC} (φ φ2 : A ⟶ B) (hca : IsCoAngular wP φ) (hst : IsPreStep wP φ)
    (hca2 : IsCoAngular wP φ2) (hst2 : IsPreStep wP φ2) (hb : wP.Base φ = wP.Base φ2)
    (α : End A) (hα : α ∈ OTri wP A) (β : End B) (hβ : β ∈ OTri wP B)
    (h : (φ ≫ (β : B ⟶ B)) = (α : A ⟶ A) ≫ φ) :
    (φ2 ≫ (β : B ⟶ B)) = (α : A ⟶ A) ≫ φ2 := by
  have hφn : wn φ = 1 := hst.1
  have hφn2 : wn φ2 = 1 := hst2.1
  have hαn : wn (α : A ⟶ A) = 1 := hα.2
  have hβn : wn (β : B ⟶ B) = 1 := hβ.2
  have hd : wd (φ ≫ (β : B ⟶ B)) = wd ((α : A ⟶ A) ≫ φ) := congrArg wd h
  simp only [wd_comp, hφn, hβn, hαn, PNat.one_coe, one_smul] at hd
  refine wHom_ext ?_ ?_
  · simp only [wd_comp, hφn2, hβn, hαn, PNat.one_coe, one_smul]
    omega
  · simp only [wn_comp, hφn2, hβn, hαn, one_mul, mul_one]

/-- **(vi)** 単元を除く忠実性。

★このモデルでは `𝒪^×(B) = {1}` なので、条件は「`φ = ψ`」に落ちる。
base-equivalent(`Vee` では自動)と metrically equivalent(`wd` が一致)、
そして両方 pre-step(`wn = 1`)から、3成分すべてが一致する。 -/
theorem wFaithfulUpToUnits {A B : wC} (φ ψ : A ⟶ B) (hbe : BaseEquivalent wP φ ψ)
    (hme : MetricallyEquivalent wP φ ψ) (hcaφ : IsCoAngular wP φ) (hstφ : IsPreStep wP φ)
    (hcaψ : IsCoAngular wP ψ) (hstψ : IsPreStep wP ψ) :
    ∃ α : End B, α ∈ OTimes wP B ∧ φ = ψ ≫ (α : B ⟶ B) := by
  refine ⟨1, (OTimes wP B).one_mem, ?_⟩
  have hd : wd φ = wd ψ := hme
  have h1 : wn φ = 1 := hstφ.1
  have h2 : wn ψ = 1 := hstψ.1
  refine wHom_ext ?_ ?_
  · show wd φ = wd (ψ ≫ (1 : End B))
    simp [hd]
  · show wn φ = wn (ψ ≫ (1 : End B))
    simp [h1, h2]

/-- **(iv)(a)** 分解の一意性。

`γ` は Frobenius 型(`wd = 0`)、`β` は pre-step(`wn = 1`、底が同型)、
`α` は pull-back(`wd = 0`、`wn = 1`)。合成は `(·, wd β, wn γ)` なので、
2つの分解から `wd β = wd β'` と `wn γ = wn γ'` が出る。
同型 `ε`、`δ` は底の同型を繋いで作る。 -/
theorem wArbFactorUniq {A B : wC} (X Y X' Y' : wC)
    (γ : A ⟶ X) (β : X ⟶ Y) (α : Y ⟶ B) (γ' : A ⟶ X') (β' : X' ⟶ Y') (α' : Y' ⟶ B)
    (hcomp : γ ≫ β ≫ α = γ' ≫ β' ≫ α')
    (hγ : IsFrobeniusType wP γ) (hβ : IsPreStep wP β) (hα : IsPullBack wP α)
    (hγ' : IsFrobeniusType wP γ') (hβ' : IsPreStep wP β') (hα' : IsPullBack wP α') :
    ∃ (δ : Y ≅ Y') (ε : X ≅ X'),
      α' = δ.inv ≫ α ∧ β' = ε.inv ≫ β ≫ δ.hom ∧ γ' = γ ≫ ε.hom := by
  haveI hbγ : IsIso (wP.Base γ) := hγ.2
  haveI hbγ' : IsIso (wP.Base γ') := hγ'.2
  haveI hbβ : IsIso (wP.Base β) := hβ.2
  haveI hbβ' : IsIso (wP.Base β') := hβ'.2
  have hdγ : wd γ = 0 := hγ.1.2
  have hdγ' : wd γ' = 0 := hγ'.1.2
  have hnβ : wn β = 1 := hβ.1
  have hnβ' : wn β' = 1 := hβ'.1
  obtain ⟨hdα, hnα⟩ := wPullBack_wd_wn α hα
  obtain ⟨hdα', hnα'⟩ := wPullBack_wd_wn α' hα'
  -- 合成の `wd` / `wn` を比べる
  have hd : wd β = wd β' := by
    have h := congrArg wd hcomp
    simp only [wd_comp, hdα, hdα', hnα, hnα', hnβ, hnβ', hdγ, hdγ',
      PNat.one_coe, one_smul, zero_add, smul_zero, add_zero] at h
    exact h
  have hn : wn γ = wn γ' := by
    have h := congrArg wn hcomp
    simp only [wn_comp, hnα, hnα', hnβ, hnβ', one_mul, mul_one] at h
    exact h
  -- 底の同型を繋ぐ
  -- 底の同型を `Iso` として繋ぐ(`inv` の合成はインスタンス解決が届かない)
  let eX : X.base ≅ X'.base :=
    (@asIso _ _ _ _ (wP.Base γ) hbγ).symm ≪≫ (@asIso _ _ _ _ (wP.Base γ') hbγ')
  let eY : Y.base ≅ Y'.base :=
    (@asIso _ _ _ _ (wP.Base β) hbβ).symm ≪≫ eX ≪≫ (@asIso _ _ _ _ (wP.Base β') hbβ')
  refine ⟨⟨⟨eY.hom, (0 : ℕ), 1⟩, ⟨eY.inv, (0 : ℕ), 1⟩,
      wHom_ext (by simp) (by simp), wHom_ext (by simp) (by simp)⟩,
    ⟨⟨eX.hom, (0 : ℕ), 1⟩, ⟨eX.inv, (0 : ℕ), 1⟩,
      wHom_ext (by simp) (by simp), wHom_ext (by simp) (by simp)⟩, ?_, ?_, ?_⟩
  · exact wHom_ext (by simp [hdα, hdα', hnα]) (by simp [hnα, hnα'])
  · exact wHom_ext (by simp [hnβ, hnβ', hd]) (by simp [hnβ, hnβ'])
  · exact wHom_ext (by simp [hdγ, hdγ']) (by simp [hn])

/-! ### ★第2段 —— 圏同値3本

★`Functor.IsEquivalence` は mathlib では **`faithful` / `full` / `essSurj` の3フィールド構造**
そのものである(2026-08-15 実測)。専用の構成補題は無く、3つを揃えれば `⟨⟩` で作れる。 -/

instance wCoaPreMul : MorphismProperty.IsMultiplicative (coaPreProp wP) :=
  coaPreProp_isMultiplicative wP wCoAngularComp

/-- `𝒞^pl-bk` の射はすべて `wd = 0`、`wn = 1` なので、hom は subsingleton。 -/
instance wPlBkSubsingleton (Z W : PlBk wP) : Subsingleton (Z ⟶ W) := by
  refine ⟨fun f g => ?_⟩
  refine InducedWideCategory.Hom.ext ?_
  exact wHom_ext (by rw [(wPullBack_wd_wn _ f.property).1, (wPullBack_wd_wn _ g.property).1])
    (by rw [(wPullBack_wd_wn _ f.property).2, (wPullBack_wd_wn _ g.property).2])

/-- **(i)(c)** `(𝒞^pl-bk)_A → 𝒟_{A_𝒟}` は圏同値。 -/
instance wPlBkFaithful (A : wC) : (plBkOverFunctor wP A).Faithful where
  map_injective {Z W} f g _ := by
    apply Over.OverMorphism.ext
    exact Subsingleton.elim _ _

instance wPlBkFull (A : wC) : (plBkOverFunctor wP A).Full where
  map_surjective {Z W} h := by
    refine ⟨Over.homMk ⟨⟨h.left, (0 : ℕ), 1⟩, wIsPullBack_of _ rfl rfl⟩ ?_, ?_⟩
    · refine InducedWideCategory.Hom.ext ?_
      refine wHom_ext ?_ ?_
      · show wd ((⟨h.left, (0 : ℕ), 1⟩ : Z.left.obj ⟶ W.left.obj) ≫ W.hom.hom) = _
        simp [(wPullBack_wd_wn _ W.hom.property).1, (wPullBack_wd_wn _ Z.hom.property).1]
      · show wn ((⟨h.left, (0 : ℕ), 1⟩ : Z.left.obj ⟶ W.left.obj) ≫ W.hom.hom) = _
        simp [(wPullBack_wd_wn _ W.hom.property).2, (wPullBack_wd_wn _ Z.hom.property).2]
    · apply Over.OverMorphism.ext
      exact Subsingleton.elim _ _

instance wPlBkEssSurj (A : wC) : (plBkOverFunctor wP A).EssSurj where
  mem_essImage Y := by
    refine ⟨Over.mk ((⟨⟨Y.hom, (0 : ℕ), 1⟩, wIsPullBack_of _ rfl rfl⟩ :
      (⟨⟨Y.left⟩⟩ : PlBk wP) ⟶ ⟨A⟩)), ⟨Over.isoMk (Iso.refl _) (Subsingleton.elim _ _)⟩⟩

theorem wPlBkEquiv (A : wC) : (plBkOverFunctor wP A).IsEquivalence :=
  ⟨inferInstance, inferInstance, inferInstance⟩

/-! #### (iii)(d) の2本 -/

/-- `𝒞^coa-pre` の射の底は同型(pre-step だから)。 -/
theorem wCoaPreBaseIso {Z W : WideSubcategory (coaPreProp wP)} (f : Z ⟶ W) :
    IsIso (wP.Base f.hom) := f.property.2.2

instance wCoaPreUnderFaithful (A : wC) : (coaPreUnderFunctor wP A).Faithful where
  map_injective {Z W} f g _ := by
    apply Under.UnderMorphism.ext
    refine InducedWideCategory.Hom.ext ?_
    haveI : Epi Z.hom.hom := wP.totEpiC _ _ _
    refine (cancel_epi Z.hom.hom).mp ?_
    have h1 : Z.hom ≫ f.right = W.hom := Under.w f
    have h2 : Z.hom ≫ g.right = W.hom := Under.w g
    calc Z.hom.hom ≫ f.right.hom = (Z.hom ≫ f.right).hom := rfl
      _ = W.hom.hom := by rw [h1]
      _ = (Z.hom ≫ g.right).hom := by rw [h2]
      _ = Z.hom.hom ≫ g.right.hom := rfl

instance wCoaPreUnderFull (A : wC) : (coaPreUnderFunctor wP A).Full where
  map_surjective {Z W} h := by
    obtain ⟨c, hc⟩ : MLe (wd Z.hom.hom) (wd W.hom.hom) := leOfHom h
    haveI hZ : IsIso (wP.Base Z.hom.hom) := Z.hom.property.2.2
    haveI hW : IsIso (wP.Base W.hom.hom) := W.hom.property.2.2
    let e : Z.right.obj.base ≅ W.right.obj.base :=
      (@asIso _ _ _ _ (wP.Base Z.hom.hom) hZ).symm ≪≫ (@asIso _ _ _ _ (wP.Base W.hom.hom) hW)
    refine ⟨Under.homMk ⟨⟨e.hom, c, 1⟩, ⟨wIsCoAngular _, rfl, (inferInstance : IsIso e.hom)⟩⟩ ?_, ?_⟩
    · refine InducedWideCategory.Hom.ext ?_
      refine wHom_ext ?_ ?_
      · show wd (Z.hom.hom ≫ (⟨e.hom, c, 1⟩ : Z.right.obj ⟶ W.right.obj)) = _
        simp only [wd_comp, wd_mk, wn_mk, PNat.one_coe, one_smul]
        omega
      · show wn (Z.hom.hom ≫ (⟨e.hom, c, 1⟩ : Z.right.obj ⟶ W.right.obj)) = _
        have h1 : wn Z.hom.hom = 1 := Z.hom.property.2.1
        have h2 : wn W.hom.hom = 1 := W.hom.property.2.1
        simp [h1, h2]
    · exact Subsingleton.elim _ _

instance wCoaPreUnderEssSurj (A : wC) : (coaPreUnderFunctor wP A).EssSurj where
  mem_essImage Y := by
    refine ⟨Under.mk ((⟨⟨𝟙 A.base, (Y : ℕ), 1⟩, ⟨wIsCoAngular _, rfl, ?_⟩⟩ :
      (⟨A⟩ : WideSubcategory (coaPreProp wP)) ⟶ ⟨A⟩)), ⟨Iso.refl _⟩⟩
    show IsIso (𝟙 A.base); infer_instance

/-- **(iii)(d)** `_A(𝒞^coa-pre) → Order(Φ(A))` は圏同値。 -/
theorem wCoaPreUnderEquiv (A : wC) : (coaPreUnderFunctor wP A).IsEquivalence :=
  ⟨inferInstance, inferInstance, inferInstance⟩

instance wCoaPreOverFaithful (A : wC) : (coaPreOverFunctor wP A).Faithful where
  map_injective {Z W} f g _ := by
    apply Over.OverMorphism.ext
    refine InducedWideCategory.Hom.ext ?_
    haveI : Mono W.hom.hom := wMono W.hom.hom
    refine (cancel_mono W.hom.hom).mp ?_
    have h1 : f.left ≫ W.hom = Z.hom := Over.w f
    have h2 : g.left ≫ W.hom = Z.hom := Over.w g
    calc f.left.hom ≫ W.hom.hom = (f.left ≫ W.hom).hom := rfl
      _ = Z.hom.hom := by rw [h1]
      _ = (g.left ≫ W.hom).hom := by rw [h2]
      _ = g.left.hom ≫ W.hom.hom := rfl

instance wCoaPreOverFull (A : wC) : (coaPreOverFunctor wP A).Full where
  map_surjective {Z W} h := by
    haveI hZ : IsIso (wP.Base Z.hom.hom) := Z.hom.property.2.2
    haveI hW : IsIso (wP.Base W.hom.hom) := W.hom.property.2.2
    obtain ⟨c, hc⟩ : MLe (wd W.hom.hom) (wd Z.hom.hom) := leOfHom h.unop
    let e : Z.left.obj.base ≅ W.left.obj.base :=
      (@asIso _ _ _ _ (wP.Base Z.hom.hom) hZ) ≪≫ (@asIso _ _ _ _ (wP.Base W.hom.hom) hW).symm
    refine ⟨Over.homMk ⟨⟨e.hom, c, 1⟩, ⟨wIsCoAngular _, rfl, (inferInstance : IsIso e.hom)⟩⟩ ?_, ?_⟩
    · refine InducedWideCategory.Hom.ext ?_
      refine wHom_ext ?_ ?_
      · show wd ((⟨e.hom, c, 1⟩ : Z.left.obj ⟶ W.left.obj) ≫ W.hom.hom) = _
        have h2 : wn W.hom.hom = 1 := W.hom.property.2.1
        simp only [wd_comp, wd_mk, h2, PNat.one_coe, one_smul]
        omega
      · show wn ((⟨e.hom, c, 1⟩ : Z.left.obj ⟶ W.left.obj) ≫ W.hom.hom) = _
        have h1 : wn Z.hom.hom = 1 := Z.hom.property.2.1
        have h2 : wn W.hom.hom = 1 := W.hom.property.2.1
        simp [h1, h2]
    · exact Subsingleton.elim _ _

instance wCoaPreOverEssSurj (A : wC) : (coaPreOverFunctor wP A).EssSurj where
  mem_essImage Y := by
    refine ⟨Over.mk ((⟨⟨𝟙 A.base, (Y.unop : ℕ), 1⟩, ⟨wIsCoAngular _, rfl, ?_⟩⟩ :
      (⟨A⟩ : WideSubcategory (coaPreProp wP)) ⟶ ⟨A⟩)), ⟨Iso.refl _⟩⟩
    show IsIso (𝟙 A.base); infer_instance

/-- **(iii)(d)** `(𝒞^coa-pre)_A → Order(Φ(A))^opp` は圏同値。 -/
theorem wCoaPreOverEquiv (A : wC) : (coaPreOverFunctor wP A).IsEquivalence :=
  ⟨inferInstance, inferInstance, inferInstance⟩

/-! ### ★自己監査で追加した2条(2026-08-15)

`(v)(b)(c)` の**一意性節**を写し落としていたので `Frobenioid.lean` に足した。
ここでその2条を `wP` について証明する。

★構造は (iv)(a) と同じ: `β` は co-angular pre-step(`wn = 1`、底が同型)、
`α` は isometric pre-step(`wd = 0`、`wn = 1`、底が同型)なので、
合成の `wd` は `wd β`、`wn` は 1。2つの分解から `wd β = wd β'` が出る。 -/

/-- 底の同型2本から `wC` の同型を作る(共通の補助)。 -/
private noncomputable def wIsoOfBases {X X' A : wC} (bX : A.base ⟶ X.base) (hX : IsIso bX)
    (bX' : A.base ⟶ X'.base) (hX' : IsIso bX') : X ≅ X' :=
  let e : X.base ≅ X'.base :=
    (@asIso _ _ _ _ bX hX).symm ≪≫ (@asIso _ _ _ _ bX' hX')
  { hom := ⟨e.hom, (0 : ℕ), 1⟩
    inv := ⟨e.inv, (0 : ℕ), 1⟩
    hom_inv_id := wHom_ext (by simp) (by simp)
    inv_hom_id := wHom_ext (by simp) (by simp) }

@[simp] private theorem wd_wIsoOfBases_hom {X X' A : wC} (bX : A.base ⟶ X.base) (hX : IsIso bX)
    (bX' : A.base ⟶ X'.base) (hX' : IsIso bX') : wd (wIsoOfBases bX hX bX' hX').hom = 0 := rfl
@[simp] private theorem wn_wIsoOfBases_hom {X X' A : wC} (bX : A.base ⟶ X.base) (hX : IsIso bX)
    (bX' : A.base ⟶ X'.base) (hX' : IsIso bX') : wn (wIsoOfBases bX hX bX' hX').hom = 1 := rfl
@[simp] private theorem wd_wIsoOfBases_inv {X X' A : wC} (bX : A.base ⟶ X.base) (hX : IsIso bX)
    (bX' : A.base ⟶ X'.base) (hX' : IsIso bX') : wd (wIsoOfBases bX hX bX' hX').inv = 0 := rfl
@[simp] private theorem wn_wIsoOfBases_inv {X X' A : wC} (bX : A.base ⟶ X.base) (hX : IsIso bX)
    (bX' : A.base ⟶ X'.base) (hX' : IsIso bX') : wn (wIsoOfBases bX hX bX' hX').inv = 1 := rfl

/-- **(v)(b)** 分解の一意性。 -/
theorem wPreStepFactorUniq {A B : wC} (X X' : wC) (β : A ⟶ X) (α : X ⟶ B)
    (β' : A ⟶ X') (α' : X' ⟶ B) (hcomp : β ≫ α = β' ≫ α')
    (hcβ : IsCoAngular wP β) (hsβ : IsPreStep wP β)
    (hiα : IsIsometric wP α) (hsα : IsPreStep wP α)
    (hcβ' : IsCoAngular wP β') (hsβ' : IsPreStep wP β')
    (hiα' : IsIsometric wP α') (hsα' : IsPreStep wP α') :
    ∃ γ : X ≅ X', α' = γ.inv ≫ α ∧ β' = β ≫ γ.hom := by
  have hnβ : wn β = 1 := hsβ.1
  have hnβ' : wn β' = 1 := hsβ'.1
  have hnα : wn α = 1 := hsα.1
  have hnα' : wn α' = 1 := hsα'.1
  have hdα : wd α = 0 := hiα
  have hdα' : wd α' = 0 := hiα'
  have hd : wd β = wd β' := by
    have h := congrArg wd hcomp
    simp only [wd_comp, hdα, hdα', hnα, hnα', PNat.one_coe, one_smul, zero_add] at h
    exact h
  refine ⟨wIsoOfBases (wP.Base β) hsβ.2 (wP.Base β') hsβ'.2, ?_, ?_⟩
  · exact wHom_ext (by simp [hdα, hdα', hnα]) (by simp [hnα, hnα'])
  · exact wHom_ext (by simp [hd]) (by simp [hnβ, hnβ'])

/-- **(v)(c)** 分解の一意性。 -/
theorem wPreStepFactorUniq' {A B : wC} (X X' : wC) (β : A ⟶ X) (α : X ⟶ B)
    (β' : A ⟶ X') (α' : X' ⟶ B) (hcomp : β ≫ α = β' ≫ α')
    (hiβ : IsIsometric wP β) (hsβ : IsPreStep wP β)
    (hcα : IsCoAngular wP α) (hsα : IsPreStep wP α)
    (hiβ' : IsIsometric wP β') (hsβ' : IsPreStep wP β')
    (hcα' : IsCoAngular wP α') (hsα' : IsPreStep wP α') :
    ∃ γ : X ≅ X', α' = γ.inv ≫ α ∧ β' = β ≫ γ.hom := by
  have hnβ : wn β = 1 := hsβ.1
  have hnβ' : wn β' = 1 := hsβ'.1
  have hnα : wn α = 1 := hsα.1
  have hnα' : wn α' = 1 := hsα'.1
  have hdβ : wd β = 0 := hiβ
  have hdβ' : wd β' = 0 := hiβ'
  have hd : wd α = wd α' := by
    have h := congrArg wd hcomp
    simp only [wd_comp, hdβ, hdβ', hnα, hnα', PNat.one_coe, smul_zero, add_zero] at h
    exact h
  refine ⟨wIsoOfBases (wP.Base β) hsβ.2 (wP.Base β') hsβ'.2, ?_, ?_⟩
  · exact wHom_ext (by simp [hd]) (by simp [hnα, hnα'])
  · exact wHom_ext (by simp [hdβ, hdβ']) (by simp [hnβ, hnβ'])


/-! ### ★★判定 —— `wP` は Frobenioid である -/

/-- ★★**`wP` は [FrdI] Definition 1.3 の意味で Frobenioid である。**

したがって **`Frobenioid` は空ではない** ——
`Definition 1.3` は退化した(誰も満たさない)条件の束ではない。

★構成要素はすべて上で個別に証明したもので、`sorry` も「手で確かめた」も含まない。 -/
theorem wIsFrobenioid : Frobenioid wP where
  core :=
    { baseSurj := wBaseSurj
      preStepSpan := wPreStepSpan
      plBkEquiv := wPlBkEquiv
      frobDegSurj := wFrobDegSurj
      frobDegUniq := wFrobDegUniq
      coAngularComp := wCoAngularComp
      coAngularOfPreStep := wCoAngularOfPreStep
      otriFwd := wOtriFwd
      otriBwd := wOtriBwd
      otriBase := wOtriBase
      arbFactor := wArbFactor
      arbFactorUniq := wArbFactorUniq
      pullBackLB := wPullBackLB
      preStepMono := fun φ _ => wMono φ
      preStepFactor := wPreStepFactor
      preStepFactorUniq := wPreStepFactorUniq
      preStepFactor' := wPreStepFactor'
      preStepFactorUniq' := wPreStepFactorUniq'
      faithfulUpToUnits := wFaithfulUpToUnits
      isotropicHullExists := wIsotropicHullExists
      isotropicClosed := fun φ h => wIsotropicClosed φ h }
  coaPreUnderEquiv := wCoaPreUnderEquiv
  coaPreOverEquiv := wCoaPreOverEquiv

/-- ★**Frobenioid は空でない**(上の系)。 -/
theorem frobenioid_nonempty : ∃ (P : PreFrobenioid wC wΦ), Frobenioid P :=
  ⟨wP, wIsFrobenioid⟩

end ABC3.Found.FrdI
