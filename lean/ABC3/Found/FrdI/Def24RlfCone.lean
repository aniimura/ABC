/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Def24Rlf
import Mathlib.RingTheory.Flat.TorsionFree
import Mathlib.RingTheory.Flat.Domain

/-!
# 実化 `Φ^rlf` —— **`ℝ ⊗_ℤ Φ^gp` の中の `ℝ≥0`-錐**として

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.51 / p.103。

原文 (FrdI p.103):
> [i.e., the “realification” of Definition 2.4, (i)] and the rational function monoid

原文 (FrdI p.103):
> Proposition 5.3. (Realifications of Frobenioids) Suppose that Φ is perf-

## ★★なぜ模型を取り替えたか

`Def24Rlf.lean` は係数拡大を **`S ⊗_ℕ M`** と書いた。`S = ℚ≥0` では
`Def24PfT.lean` が条件 (a) を証明したが、`S = ℝ≥0` では
**半環上の平坦性**(`ℝ≥0` が `ℚ≥0` 上平坦)が要り、これは
「`ℝ≥0` が有限生成自由 `ℚ≥0`-半加群の濾過余極限である」という
凸幾何の議論を通す必要がある(依存グラフの鎖 `rlf` の節点 `rlf-flat`)。

★★★**そこを迂回する** —— **群化してから実係数にする**:

  `Φ^rlf(A) := ` (`ℝ ⊗_ℤ Φ(A)^gp` の中で `Φ(A)` の像が生成する `ℝ≥0`-錐)

こうすると条件 (a) の**第 1 成分が無料になる**:

* `Φ(A)^gp → Φ(B)^gp` は単射(`gpMap_injective`、`Φ` が integral)
* **`ℝ` は `ℤ` 上平坦**(mathlib の instance で出る。`ℤ` は PID で `ℝ` は捻れ無し)
* ゆえに `ℝ ⊗_ℤ Φ(A)^gp → ℝ ⊗_ℤ Φ(B)^gp` も単射
* 錐への制限も単射

★★半環上の平坦性が要らなくなり、**環上の平坦性(mathlib 在庫)**で済む。

## ★★逸脱(記録)

★これは `Def24Rlf.lean` の `ScT ℝ≥0 M = ℝ≥0 ⊗_ℕ M` とは**別の模型**である
(自然な全射 `ℝ≥0 ⊗_ℕ M ↠ (錐)` があるが、単射性はまさに `rlf-flat` である)。
★★**perf-factorial な `M` では両者は一致する**はずである ——
`Definition 2.4, (i)` の素点分解 `M^rlf ≅ ⊕_p ℝ≥0` と、
`Gp M ⊆ ⊕_p ℚ` から `ℝ ⊗_ℤ Gp M ≅ ⊕_p ℝ` で錐が `⊕_p ℝ≥0` になるからである。
★**この一致はまだ証明していない**(鎖 `rlf` の節点 `rlf-agree`)。
`Proposition 5.3` は `Φ` が perf-factorial であることを**仮定している**ので、
そこでは両模型のどちらを採っても主張は同じはずである。

## ★残る条件

★条件 (a) の**第 2 成分**(`charMap` の単射性)は、錐が **sharp** であれば
`charMap_injective_of_sharp` で出る。★錐の sharp 性は
`ModelData.Hyp` の `divisorial` が**どのみち要求する**ものなので、
新たな負担ではない。
★これは「`Σ c_k x_k = 0`(`c_k > 0` 実数、`x_k ∈ Gp M`)なら
正の**整数**係数の関係式がある」という有理性の言明に帰着する
(有限生成部分群を取れば整数行列の核の有理性)。
-/

namespace ABC3.Found.FrdI

open CategoryTheory
open scoped TensorProduct

universe v u w

variable {M N : Type w} [AddCommMonoid M] [AddCommMonoid N]

/-! ## ★1. 台 —— `ℝ ⊗_ℤ M^gp` -/

/-- ★★**実化の台** `ℝ ⊗_ℤ M^gp`。 -/
abbrev RlfV (M : Type w) [AddCommMonoid M] : Type w := ℝ ⊗[ℤ] (Gp M)

theorem toGp_zero : toGp M 0 = 0 := AddLocalization.mk_zero

theorem toGp_add (m n : M) : toGp M (m + n) = toGp M m + toGp M n := by
  show AddLocalization.mk (m + n) 0 = AddLocalization.mk m 0 + AddLocalization.mk n 0
  rw [AddLocalization.mk_add]
  congr 1
  exact (add_zero _).symm

/-- ★`M → ℝ ⊗_ℤ M^gp`(`m ↦ 1 ⊗ [m]`)。 -/
noncomputable def toRlfV : M →+ RlfV M where
  toFun m := (1 : ℝ) ⊗ₜ[ℤ] toGp M m
  map_zero' := by
    show (1 : ℝ) ⊗ₜ[ℤ] toGp M 0 = 0
    rw [toGp_zero, TensorProduct.tmul_zero]
  map_add' m n := by
    show (1 : ℝ) ⊗ₜ[ℤ] toGp M (m + n)
      = (1 : ℝ) ⊗ₜ[ℤ] toGp M m + (1 : ℝ) ⊗ₜ[ℤ] toGp M n
    rw [toGp_add, TensorProduct.tmul_add]

/-- ★★**関手性** —— `ℝ ⊗_ℤ (−)^gp`。 -/
noncomputable def rlfVMap (f : M →+ N) : RlfV M →ₗ[ℤ] RlfV N :=
  LinearMap.lTensor ℝ (gpMap M f).toIntLinearMap

@[simp] theorem rlfVMap_tmul (f : M →+ N) (r : ℝ) (g : Gp M) :
    rlfVMap f (r ⊗ₜ[ℤ] g) = r ⊗ₜ[ℤ] gpMap M f g := rfl

/-- ★★★★★**条件 (a) の第 1 成分が無料になる** —— `ℝ` は `ℤ` 上平坦。 -/
theorem rlfVMap_injective [IsCancelAdd M] [IsCancelAdd N] {f : M →+ N}
    (hf : Function.Injective f) : Function.Injective (rlfVMap f) :=
  Module.Flat.lTensor_preserves_injective_linearMap (M := ℝ)
    (gpMap M f).toIntLinearMap (gpMap_injective M hf)

/-! ## ★2. `ℝ≥0`-錐 -/

/-- ★`M` の像の非負実数倍からなる集合(錐の生成集合)。 -/
def rlfGen (M : Type w) [AddCommMonoid M] : Set (RlfV M) :=
  {x | ∃ (r : ℝ), 0 ≤ r ∧ ∃ m : M, x = r ⊗ₜ[ℤ] toGp M m}

/-- ★★**実化** `M^rlf` —— `ℝ ⊗_ℤ M^gp` の中で `M` の像が生成する `ℝ≥0`-錐。 -/
def rlfCone (M : Type w) [AddCommMonoid M] : AddSubmonoid (RlfV M) :=
  AddSubmonoid.closure (rlfGen M)

theorem toRlfV_mem_rlfCone (m : M) : toRlfV m ∈ rlfCone M :=
  AddSubmonoid.subset_closure ⟨1, zero_le_one, m, rfl⟩

theorem rlfVMap_rlfGen (f : M →+ N) {x : RlfV M} (hx : x ∈ rlfGen M) :
    rlfVMap f x ∈ rlfGen N := by
  obtain ⟨r, hr, m, rfl⟩ := hx
  exact ⟨r, hr, f m, by rw [rlfVMap_tmul, gpMap_toGp]⟩

theorem rlfVMap_mem_rlfCone (f : M →+ N) {x : RlfV M} (hx : x ∈ rlfCone M) :
    rlfVMap f x ∈ rlfCone N := by
  have hle : rlfCone M ≤ AddSubmonoid.comap ((rlfVMap f).toAddMonoidHom) (rlfCone N) :=
    AddSubmonoid.closure_le.mpr
      (fun y hy => AddSubmonoid.subset_closure (rlfVMap_rlfGen f hy))
  exact hle hx

/-- ★★**錐の間の写像**。 -/
noncomputable def rlfConeMap (f : M →+ N) : rlfCone M →+ rlfCone N :=
  AddMonoidHom.codRestrict
    (((rlfVMap f).toAddMonoidHom).comp (rlfCone M).subtype) _
    (fun x => rlfVMap_mem_rlfCone f x.2)

@[simp] theorem rlfConeMap_coe (f : M →+ N) (x : rlfCone M) :
    ((rlfConeMap f x : rlfCone N) : RlfV N) = rlfVMap f (x : RlfV M) := rfl

theorem rlfConeMap_injective [IsCancelAdd M] [IsCancelAdd N] {f : M →+ N}
    (hf : Function.Injective f) : Function.Injective (rlfConeMap f) := by
  intro x y h
  exact Subtype.ext (rlfVMap_injective hf
    (congrArg (fun t : rlfCone N => (t : RlfV N)) h))

theorem rlfConeMap_surjective {f : M →+ N} (hf : Function.Surjective f) :
    Function.Surjective (rlfConeMap f) := by
  intro y
  have hle : rlfCone N ≤ AddSubmonoid.map ((rlfVMap f).toAddMonoidHom) (rlfCone M) := by
    refine AddSubmonoid.closure_le.mpr (fun z hz => ?_)
    obtain ⟨r, hr, n, rfl⟩ := hz
    obtain ⟨m, rfl⟩ := hf n
    refine ⟨r ⊗ₜ[ℤ] toGp M m, AddSubmonoid.subset_closure ⟨r, hr, m, rfl⟩, ?_⟩
    show rlfVMap f (r ⊗ₜ[ℤ] toGp M m) = r ⊗ₜ[ℤ] toGp N (f m)
    rw [rlfVMap_tmul, gpMap_toGp]
  obtain ⟨x, hx, hxy⟩ := hle y.2
  exact ⟨⟨x, hx⟩, Subtype.ext hxy⟩

/-! ## ★3. `Φ^rlf : MonoidOn 𝒟` -/

theorem rlfVMap_id (x : RlfV M) : rlfVMap (AddMonoidHom.id M) x = x := by
  have h : gpMap M (AddMonoidHom.id M) = AddMonoidHom.id (Gp M) := gpMap_id M
  show LinearMap.lTensor ℝ (gpMap M (AddMonoidHom.id M)).toIntLinearMap x = x
  rw [h]
  show LinearMap.lTensor ℝ (LinearMap.id : Gp M →ₗ[ℤ] Gp M) x = x
  rw [LinearMap.lTensor_id]
  rfl

theorem rlfVMap_comp {O : Type w} [AddCommMonoid O] (f : M →+ N) (g : N →+ O) (x : RlfV M) :
    rlfVMap (g.comp f) x = rlfVMap g (rlfVMap f x) := by
  have h : gpMap M (g.comp f) = (gpMap N g).comp (gpMap M f) := gpMap_comp M f g
  show LinearMap.lTensor ℝ (gpMap M (g.comp f)).toIntLinearMap x = _
  rw [h]
  show LinearMap.lTensor ℝ
    ((gpMap N g).toIntLinearMap ∘ₗ (gpMap M f).toIntLinearMap) x = _
  rw [LinearMap.lTensor_comp]
  rfl

theorem rlfConeMap_id (x : rlfCone M) : rlfConeMap (AddMonoidHom.id M) x = x :=
  Subtype.ext (rlfVMap_id (x : RlfV M))

theorem rlfConeMap_comp {O : Type w} [AddCommMonoid O] (f : M →+ N) (g : N →+ O)
    (x : rlfCone M) :
    rlfConeMap (g.comp f) x = rlfConeMap g (rlfConeMap f x) :=
  Subtype.ext (rlfVMap_comp f g (x : RlfV M))

/-! ## ★★3. 錐が sharp であるための十分条件

★★錐の sharp 性は `ModelData.Hyp` の `divisorial` がどのみち要求するものである。
★**`M` 上で狭義正の線型汎関数があれば十分**であり、
幾何（`Example 6.1`）では係数の総和、
算術（`Example 6.3`）では次数がそれである。 -/

/-- ★`ℓ : M^gp →+ ℝ` を `ℝ ⊗_ℤ M^gp →ₗ[ℤ] ℝ` へ延ばす。 -/
noncomputable def rlfLift (l : Gp M →+ ℝ) : RlfV M →ₗ[ℤ] ℝ :=
  TensorProduct.lift
    { toFun := fun r => r • l.toIntLinearMap
      map_add' := by intro a b; ext g; simp [add_smul]
      map_smul' := by intro n a; ext g; simp [mul_assoc] }

@[simp] theorem rlfLift_tmul (l : Gp M →+ ℝ) (r : ℝ) (g : Gp M) :
    rlfLift l (r ⊗ₜ[ℤ] g) = r * l g := by
  show (r • l.toIntLinearMap) g = r * l g
  simp [smul_eq_mul]

/-- ★★★★★**錐が sharp であるための十分条件** —— `M` 上で狭義正の線型汎関数。

★証明は `AddSubmonoid.closure_induction` で「`0 ≤ L x` かつ `L x = 0 → x = 0`」を
錐全体へ伸ばすだけである（この述語は和で閉じている）。 -/
theorem isSharp_rlfCone_of_pos (l : Gp M →+ ℝ)
    (hpos : ∀ m : M, m ≠ 0 → 0 < l (toGp M m)) : IsSharp (rlfCone M) := by
  set L := rlfLift l with hL
  have key : ∀ x ∈ rlfCone M, 0 ≤ L x ∧ (L x = 0 → x = 0) := by
    intro x hx
    refine AddSubmonoid.closure_induction ?_ ?_ ?_ hx
    · rintro y ⟨r, hr, m, rfl⟩
      by_cases hm : m = 0
      · subst hm
        rw [toGp_zero, TensorProduct.tmul_zero]
        simp
      · have hlm : 0 < l (toGp M m) := hpos m hm
        refine ⟨by rw [hL, rlfLift_tmul]; positivity, ?_⟩
        intro h0
        rw [hL, rlfLift_tmul] at h0
        have hr0 : r = 0 := by
          rcases mul_eq_zero.mp h0 with h | h
          · exact h
          · exact absurd h hlm.ne'
        rw [hr0, TensorProduct.zero_tmul]
    · exact ⟨by simp, fun _ => rfl⟩
    · rintro a b _ _ ⟨ha0, ha1⟩ ⟨hb0, hb1⟩
      refine ⟨by rw [map_add]; linarith, ?_⟩
      intro h
      rw [map_add] at h
      have hA : L a = 0 := by linarith
      have hB : L b = 0 := by linarith
      rw [ha1 hA, hb1 hB, add_zero]
  intro a ha
  obtain ⟨u, rfl⟩ := ha
  have h1 := key ((u.val : rlfCone M) : RlfV M) u.val.2
  have h2 := key ((u.neg : rlfCone M) : RlfV M) u.neg.2
  have hsum : ((u.val : rlfCone M) : RlfV M) + ((u.neg : rlfCone M) : RlfV M) = 0 :=
    congrArg (fun t : rlfCone M => (t : RlfV M)) u.val_neg
  have h3 : L ((u.val : rlfCone M) : RlfV M) + L ((u.neg : rlfCone M) : RlfV M) = 0 := by
    rw [← map_add, hsum, map_zero]
  have h4 : L ((u.val : rlfCone M) : RlfV M) = 0 := by linarith [h1.1, h2.1]
  exact Subtype.ext (h1.2 h4)

/-- ★★★★★★**錐が sharp であるための十分条件(族版)** ——
`M` 上で**非負**な線型汎関数の族であって、`M` の点を**分離する**もの。

★★★これが `perf-factorial` にそのまま効く —— 素点ごとの `ord_p`
(`factorMap ι · p`)は**非負**であり、`factorMap` の単射性がちょうど**分離**である。
★★1 本の狭義正な汎関数(`isSharp_rlfCone_of_pos`)を要求すると
素点が無限個のとき「係数の総和」が有限和にならないが、
**族にすれば有限和は要らない** —— 錐の生成元 `r ⊗ toGp m` は
1 点 `m` しか含まないので、各素点で別々に消えることを見ればよい。 -/
theorem isSharp_rlfCone_of_family {ι : Type*} (l : ι → (Gp M →+ ℝ))
    (hnonneg : ∀ (i : ι) (m : M), 0 ≤ l i (toGp M m))
    (hsep : ∀ m : M, (∀ i : ι, l i (toGp M m) = 0) → m = 0) :
    IsSharp (rlfCone M) := by
  have key : ∀ x ∈ rlfCone M, (∀ i, 0 ≤ rlfLift (l i) x) ∧
      ((∀ i, rlfLift (l i) x = 0) → x = 0) := by
    intro x hx
    refine AddSubmonoid.closure_induction ?_ ?_ ?_ hx
    · rintro y ⟨r, hr, m, rfl⟩
      refine ⟨fun i => ?_, ?_⟩
      · rw [rlfLift_tmul]
        exact mul_nonneg hr (hnonneg i m)
      · intro h0
        by_cases hr0 : r = 0
        · rw [hr0, TensorProduct.zero_tmul]
        · have hm : m = 0 := hsep m fun i => by
            have hi := h0 i
            rw [rlfLift_tmul] at hi
            rcases mul_eq_zero.mp hi with h | h
            · exact absurd h hr0
            · exact h
          rw [hm, toGp_zero, TensorProduct.tmul_zero]
    · exact ⟨fun _ => by simp, fun _ => rfl⟩
    · rintro a b _ _ ⟨ha0, ha1⟩ ⟨hb0, hb1⟩
      refine ⟨fun i => by rw [map_add]; linarith [ha0 i, hb0 i], ?_⟩
      intro h
      have hA : ∀ i, rlfLift (l i) a = 0 := fun i => by
        have hi := h i
        rw [map_add] at hi
        linarith [ha0 i, hb0 i]
      have hB : ∀ i, rlfLift (l i) b = 0 := fun i => by
        have hi := h i
        rw [map_add] at hi
        linarith [ha0 i, hb0 i]
      rw [ha1 hA, hb1 hB, add_zero]
  intro a ha
  obtain ⟨u, rfl⟩ := ha
  have h1 := key ((u.val : rlfCone M) : RlfV M) u.val.2
  have h2 := key ((u.neg : rlfCone M) : RlfV M) u.neg.2
  have hsum : ((u.val : rlfCone M) : RlfV M) + ((u.neg : rlfCone M) : RlfV M) = 0 :=
    congrArg (fun t : rlfCone M => (t : RlfV M)) u.val_neg
  refine Subtype.ext (h1.2 fun i => ?_)
  have h3 : rlfLift (l i) ((u.val : rlfCone M) : RlfV M)
      + rlfLift (l i) ((u.neg : rlfCone M) : RlfV M) = 0 := by
    rw [← map_add, hsum, map_zero]
  linarith [h1.1 i, h2.1 i]

variable {D : Type u} [Category.{v} D]

/-- ★`Φ^rlf` の台となる反変関手。 -/
noncomputable def rlfConeFunctor (Φ : MonoidOn.{v, u, w} D) : Dᵒᵖ ⥤ AddCommMonCat.{w} where
  obj X := AddCommMonCat.of (rlfCone (Φ.val X.unop))
  map f := AddCommMonCat.ofHom (rlfConeMap (Φ.map f.unop))
  map_id X := by
    refine AddCommMonCat.ext (fun x => ?_)
    have h : Φ.map (𝟙 X.unop) = AddMonoidHom.id _ := by ext a; exact Φ.map_id _ a
    show rlfConeMap (Φ.map (𝟙 X.unop)) x = _
    rw [h]
    exact rlfConeMap_id x
  map_comp {X Y Z} f g := by
    refine AddCommMonCat.ext (fun x => ?_)
    have h : Φ.map (g.unop ≫ f.unop) = (Φ.map g.unop).comp (Φ.map f.unop) := by
      ext a; exact Φ.map_comp _ _ a
    show rlfConeMap (Φ.map (g.unop ≫ f.unop)) x = _
    rw [h]
    exact rlfConeMap_comp _ _ x

/-- ★★★★★★**`Φ^rlf` を `𝒟` 上の単系として立てる**。

★条件 (a) の**第 1 成分は無料**(`ℝ` の `ℤ` 上の平坦性)。
★残るのは**錐が sharp** であること —— これは `ModelData.Hyp` の
`divisorial` がどのみち要求するものである。 -/
noncomputable def phiRlfConeOn (Φ : MonoidOn.{v, u, w} D)
    (hcancel : ∀ A : D, IsCancelAdd (Φ.val A))
    (hsharp : ∀ A : D, IsSharp (rlfCone (Φ.val A))) : MonoidOn.{v, u, w} D where
  functor := rlfConeFunctor Φ
  charInj {A B} α := by
    haveI := hcancel A
    haveI := hcancel B
    have hinj : Function.Injective (rlfConeMap (Φ.map α)) :=
      rlfConeMap_injective (Φ.map_injective α)
    exact ⟨hinj, charMap_injective_of_sharp (hsharp B) hinj⟩
  fsmIso {A B} α hα := by
    haveI := hcancel A
    haveI := hcancel B
    exact ⟨rlfConeMap_injective (Φ.map_injective α),
      rlfConeMap_surjective (Φ.fsmIso α hα).2⟩

/-- ★divisorial なら簡約性は自動。 -/
noncomputable def phiRlfConeOnOfDivisorial (Φ : MonoidOn.{v, u, w} D)
    (hdiv : Φ.IsDivisorialOn)
    (hsharp : ∀ A : D, IsSharp (rlfCone (Φ.val A))) : MonoidOn.{v, u, w} D :=
  phiRlfConeOn Φ (fun A => isCancelAdd_of_isIntegralMonoid _ (hdiv A).1.1) hsharp

/-- ★★★★★★★**各 `Φ(A)` に狭義正の線型汎関数があれば、
`Φ^rlf : MonoidOn D` は仮定なしで立つ**。

★幾何（`Example 6.1`）では係数の総和、算術（`Example 6.3`）では
次数がその汎関数である。 -/
noncomputable def phiRlfConeOnOfPos (Φ : MonoidOn.{v, u, w} D)
    (hdiv : Φ.IsDivisorialOn)
    (l : ∀ A : D, Gp (Φ.val A) →+ ℝ)
    (hpos : ∀ (A : D) (m : Φ.val A), m ≠ 0 → 0 < l A (toGp (Φ.val A) m)) :
    MonoidOn.{v, u, w} D :=
  phiRlfConeOnOfDivisorial Φ hdiv (fun A => isSharp_rlfCone_of_pos (l A) (hpos A))

/-! ### ★出典の紐付け -/

/-- ★locator —— `Definition 2.4, (i)` の実化(★**錐模型**。
`Def24.lean` の素点分解版・`Def24Rlf.lean` のテンソル版とは別の模型であり、
一致は未証明)。 -/
def rlfCone.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 51, item := "Definition 2.4, (i) — realification(錐模型)",
    sectionId := "frdi-def-2-4" }

/-- ★locator —— `Φ^rlf` が `𝒟` 上の単系であること。 -/
def phiRlfConeOn.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 103, item := "Proposition 5.3 — Φ^rlf は monoid on 𝒟(錐模型)",
    sectionId := "frdi-prop-5-3" }

end ABC3.Found.FrdI
