/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Def24Rlf
import ABC3.Found.FrdI.Prop32

/-!
# `M^pf ≃ ℚ≥0 ⊗_ℕ M` —— 完全化の**商模型とテンソル模型の一致**

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.19。

原文 (FrdI p.19):
> is a monoid on D, then Φ determines monoids “Φchar”, “Φgp”, Φpf” on D [i.e.,

## ★★このファイルが外す障壁

`Def24Rlf.lean` は係数拡大を `M ⊗_S := S ⊗_ℕ M` と書き、
`S = ℝ≥0` が実化・`S = ℚ≥0` が完全化になるようにした。
★そこで `MonoidOn` の条件 (a)(characteristic injectivity)は
**テンソルの平坦性**にあたり、一般の可換単系では成り立たないので
`phiScOn` は**それを仮定 `hcharInj` として受け取る**形にしてあった
(`ResearchPaper/frdi-decomposition.json` の鎖 `rlf` の節点 `rlf-flat`)。

★★★**`S = ℚ≥0` の場合、それは証明できる。**
平坦性を直接示すのではなく、**商模型 `M^pf`(`MonoidVocabulary.lean` の `Pf`)と
同一視する**のである:

* `pfToPfT : M^pf → ℚ≥0 ⊗_ℕ M`、`m/a ↦ (1/a) ⊗ m`
* `pfTToPf : ℚ≥0 ⊗_ℕ M → M^pf`、`q ⊗ m ↦ (q.num · m)/q.den`
* この 2 本は互いに逆(`pfEquivPfT`)で、**自然**である(`pfToPfT_natural`)

`Pf` の側では単射性の保存が **5 行**で出る(`Pf.map_injective` ——
`f((k a')·m) = f((k a)·m')` から単射性で `(k a')·m = (k a)·m'`)。
★したがって `ℚ≥0 ⊗_ℕ −` も単射性を保つ(`scMap_injective_nnrat`)。

条件 (a) の第 2 成分(`charMap` の単射性)は、`Φ` が divisorial なら
各 `Φ(A)` が sharp、ゆえに `M^pf` も sharp(`Pf.isSharp_pf`)、
ゆえに `ℚ≥0 ⊗_ℕ M` も sharp(`isSharp_pfT`)なので
`charMap_injective_of_sharp` で出る。

★★**結論**: `phiPfTOn`(下)は**仮定なしで** `Φ^pf : MonoidOn 𝒟` を与える。
★残る `rlf-flat` は **`S = ℝ≥0` の場合だけ**である
(`ℝ≥0` には `Pf` のような商模型が無い —— `ℝ≥0` は `ℚ≥0` 上の
有限生成自由半加群の濾過余極限であり、そこを通す必要がある)。

## ★在庫

`MonoidOn.pfOn`(`ElementaryFrobenioid.lean`)は**商模型 `Pf` による** `Φ^pf` を
すでに与えている。本ファイルはそれを**テンソル模型に移す**もので、
`Prop53Rlf.lean` の `scModel` / `CuntrPf`(係数 `S` で書かれている)が
そのまま使えるようにするためのものである。
-/

namespace ABC3.Found.FrdI

open CategoryTheory
open scoped TensorProduct

universe v u w

variable {M N : Type w} [AddCommMonoid M] [AddCommMonoid N]

/-! ## ★1. 分母まわりの小道具 -/

/-- ★`1/a ∈ ℚ≥0`。 -/
noncomputable def pnatInv (a : ℕ+) : ℚ≥0 := ((a : ℕ) : ℚ≥0)⁻¹

theorem pnatInv_ne_zero (a : ℕ+) : pnatInv a ≠ 0 := by
  rw [pnatInv, ne_eq, inv_eq_zero]
  exact_mod_cast a.ne_zero

theorem pnatInv_mul (a : ℕ+) : pnatInv a * (((a : ℕ+) : ℕ) : ℚ≥0) = ((1 : ℕ) : ℚ≥0) := by
  have ha : (((a : ℕ+) : ℕ) : ℚ≥0) ≠ 0 := by exact_mod_cast a.ne_zero
  rw [pnatInv, Nat.cast_one]
  field_simp

/-- ★`ℚ≥0` の元の分母を `ℕ+` として。 -/
def nnratDen (q : ℚ≥0) : ℕ+ := ⟨q.den, q.den_pos⟩

@[simp] theorem nnratDen_coe (q : ℚ≥0) : ((nnratDen q : ℕ+) : ℕ) = q.den := rfl

/-- ★`q.num = q · q.den`(`ℚ≥0` の中で)。 -/
theorem nnrat_num_eq (q : ℚ≥0) : ((q.num : ℕ) : ℚ≥0) = q * ((q.den : ℕ) : ℚ≥0) := by
  have hq : ((q.den : ℕ) : ℚ≥0) ≠ 0 := by exact_mod_cast q.den_pos.ne'
  have h := q.num_div_den
  field_simp at h
  rw [h, mul_comm]

/-- ★和の分子・分母の**交叉等式**。 -/
theorem nnrat_add_cross (q q' : ℚ≥0) :
    (q.den * q'.den) * (q + q').num = (q + q').den * (q'.den * q.num + q.den * q'.num) := by
  have h : ((q.den * q'.den : ℕ) : ℚ≥0) * (((q + q').num : ℕ) : ℚ≥0)
      = (((q + q').den : ℕ) : ℚ≥0) * ((q'.den * q.num + q.den * q'.num : ℕ) : ℚ≥0) := by
    push_cast
    rw [nnrat_num_eq q, nnrat_num_eq q', nnrat_num_eq (q + q')]
    ring
  exact_mod_cast h

/-- ★`ℕ` 倍の分子・分母の交叉等式。 -/
theorem nnrat_nsmul_cross (n : ℕ) (q : ℚ≥0) :
    q.den * (n • q).num = (n • q).den * (n * q.num) := by
  have h : ((q.den : ℕ) : ℚ≥0) * (((n • q).num : ℕ) : ℚ≥0)
      = (((n • q).den : ℕ) : ℚ≥0) * ((n * q.num : ℕ) : ℚ≥0) := by
    push_cast
    rw [nnrat_num_eq q, nnrat_num_eq (n • q), nsmul_eq_mul]
    ring
  exact_mod_cast h

/-! ## ★2. `Pf` の小道具 -/

theorem Pf.mk_zero (a : ℕ+) : Pf.mk (0 : M) a = 0 := by
  refine Pf.sound 1 ?_
  simp only [smul_zero]

/-- ★★分子・分母の**交叉等式**から `M^pf` の等式が出る。 -/
theorem Pf.mk_nsmul_eq_of_cross' {N P : ℕ} {D Q : ℕ+} (m : M)
    (h : (Q : ℕ) * N = (D : ℕ) * P) :
    Pf.mk (N • m) D = Pf.mk (P • m) Q := by
  refine Pf.sound 1 ?_
  simp only [PNat.one_coe, one_mul, smul_smul]
  rw [h]

/-! ## ★3. `M^pf → ℚ≥0 ⊗_ℕ M` -/

/-- ★`m/a ↦ (1/a) ⊗ m`(台の関数)。 -/
noncomputable def pfToPfTFun (x : Pf M) : PfT M :=
  Quotient.liftOn x (fun y : M × ℕ+ => pnatInv y.2 ⊗ₜ[ℕ] y.1)
    (by
      rintro ⟨m, a⟩ ⟨m', a'⟩ ⟨k, e⟩
      have key : ∀ (b c : ℕ+) (z : M), pnatInv b ⊗ₜ[ℕ] z
          = pnatInv (k * b * c) ⊗ₜ[ℕ] (((k : ℕ) * (c : ℕ)) • z) := by
        intro b c z
        rw [← TensorProduct.smul_tmul]
        congr 1
        have hb : ((b : ℕ) : ℚ≥0) ≠ 0 := by exact_mod_cast b.ne_zero
        have hk : ((k : ℕ) : ℚ≥0) ≠ 0 := by exact_mod_cast k.ne_zero
        have hc : ((c : ℕ) : ℚ≥0) ≠ 0 := by exact_mod_cast c.ne_zero
        rw [nsmul_eq_mul, pnatInv, pnatInv]
        push_cast
        field_simp
      have h1 : pnatInv a ⊗ₜ[ℕ] m = pnatInv (k * a * a') ⊗ₜ[ℕ] (((k : ℕ) * (a' : ℕ)) • m) :=
        key a a' m
      have h2 : pnatInv a' ⊗ₜ[ℕ] m' = pnatInv (k * a' * a) ⊗ₜ[ℕ] (((k : ℕ) * (a : ℕ)) • m') :=
        key a' a m'
      rw [h1, h2, e]
      congr 2
      exact (mul_right_comm k a a').symm ▸ rfl)

@[simp] theorem pfToPfTFun_mk (m : M) (a : ℕ+) :
    pfToPfTFun (Pf.mk m a) = pnatInv a ⊗ₜ[ℕ] m := rfl

/-- ★★**`M^pf → ℚ≥0 ⊗_ℕ M`** —— `m/a ↦ (1/a) ⊗ m`。 -/
noncomputable def pfToPfT : Pf M →+ PfT M where
  toFun := pfToPfTFun
  map_zero' := by
    show pfToPfTFun (Pf.mk (0 : M) 1) = 0
    rw [pfToPfTFun_mk, TensorProduct.tmul_zero]
  map_add' x y := by
    induction x using Pf.inductionOn with | _ m a =>
    induction y using Pf.inductionOn with | _ m' a' =>
    show pfToPfTFun (Pf.mk ((a' : ℕ) • m + (a : ℕ) • m') (a * a'))
      = pnatInv a ⊗ₜ[ℕ] m + pnatInv a' ⊗ₜ[ℕ] m'
    rw [pfToPfTFun_mk, TensorProduct.tmul_add, ← TensorProduct.smul_tmul,
      ← TensorProduct.smul_tmul]
    congr 2
    · have ha : ((a : ℕ) : ℚ≥0) ≠ 0 := by exact_mod_cast a.ne_zero
      have ha' : ((a' : ℕ) : ℚ≥0) ≠ 0 := by exact_mod_cast a'.ne_zero
      rw [nsmul_eq_mul, pnatInv, pnatInv]
      push_cast
      field_simp
    · have ha : ((a : ℕ) : ℚ≥0) ≠ 0 := by exact_mod_cast a.ne_zero
      have ha' : ((a' : ℕ) : ℚ≥0) ≠ 0 := by exact_mod_cast a'.ne_zero
      rw [nsmul_eq_mul, pnatInv, pnatInv]
      push_cast
      field_simp

/-! ## ★4. `ℚ≥0 ⊗_ℕ M → M^pf` -/

/-- ★`q ⊙ m := (q.num · m)/q.den`。 -/
noncomputable def qsmulPf (q : ℚ≥0) : M →+ Pf M where
  toFun m := Pf.mk (q.num • m) (nnratDen q)
  map_zero' := by rw [smul_zero, Pf.mk_zero]
  map_add' m n := by rw [smul_add, ← Pf.mk_add_mk_same]

@[simp] theorem qsmulPf_apply (q : ℚ≥0) (m : M) :
    qsmulPf q m = Pf.mk (q.num • m) (nnratDen q) := rfl

/-- ★`q ↦ (m ↦ q ⊙ m)` は加法的(和の交叉等式による)。 -/
noncomputable def qsmulPfHom : ℚ≥0 →+ (M →+ Pf M) where
  toFun := qsmulPf
  map_zero' := by
    ext m
    show Pf.mk ((0 : ℚ≥0).num • m) (nnratDen 0) = 0
    rw [show (0 : ℚ≥0).num = 0 from rfl, zero_smul, Pf.mk_zero]
  map_add' q q' := by
    ext m
    show Pf.mk ((q + q').num • m) (nnratDen (q + q'))
      = Pf.mk (q.num • m) (nnratDen q) + Pf.mk (q'.num • m) (nnratDen q')
    rw [Pf.mk_add_mk]
    have hrw : ((nnratDen q' : ℕ+) : ℕ) • (q.num • m) + ((nnratDen q : ℕ+) : ℕ) • (q'.num • m)
        = (q'.den * q.num + q.den * q'.num) • m := by
      rw [add_smul, smul_smul, smul_smul]
      rfl
    rw [hrw]
    refine Pf.mk_nsmul_eq_of_cross' m ?_
    show (q.den * q'.den) * (q + q').num
      = (q + q').den * (q'.den * q.num + q.den * q'.num)
    exact nnrat_add_cross q q'

/-- ★★**`ℚ≥0 ⊗_ℕ M → M^pf`**。 -/
noncomputable def pfTToPf : PfT M →+ Pf M :=
  TensorProduct.liftAddHom (R := ℕ) (qsmulPfHom (M := M)) (by
    intro n q m
    show Pf.mk ((n • q).num • m) (nnratDen (n • q))
      = Pf.mk (q.num • (n • m)) (nnratDen q)
    rw [smul_smul]
    refine Pf.mk_nsmul_eq_of_cross' m ?_
    show q.den * (n • q).num = (n • q).den * (q.num * n)
    rw [mul_comm q.num n]
    exact nnrat_nsmul_cross n q)

@[simp] theorem pfTToPf_tmul (q : ℚ≥0) (m : M) :
    pfTToPf (q ⊗ₜ[ℕ] m) = Pf.mk (q.num • m) (nnratDen q) := rfl

/-- ★`q = n/d` と書けるときの `pfTToPf` の値。 -/
theorem pfTToPf_tmul_of (n : ℕ) (d : ℕ+) (q : ℚ≥0)
    (hq : q * (((d : ℕ+) : ℕ) : ℚ≥0) = (n : ℚ≥0)) (m : M) :
    pfTToPf (q ⊗ₜ[ℕ] m) = Pf.mk (n • m) d := by
  rw [pfTToPf_tmul]
  refine Pf.mk_nsmul_eq_of_cross' m ?_
  show ((d : ℕ+) : ℕ) * q.num = q.den * n
  have h : (((d : ℕ+) : ℕ) : ℚ≥0) * ((q.num : ℕ) : ℚ≥0)
      = ((q.den : ℕ) : ℚ≥0) * ((n : ℕ) : ℚ≥0) := by
    rw [nnrat_num_eq q, ← hq]
    ring
  exact_mod_cast h

/-! ## ★5. 2 本は互いに逆 -/

/-- ★★**往復 1** —— `M^pf → ℚ≥0 ⊗ M → M^pf` は恒等。 -/
theorem pfTToPf_pfToPfT (x : Pf M) : pfTToPf (pfToPfT x) = x := by
  induction x using Pf.inductionOn with | _ m a =>
  show pfTToPf (pnatInv a ⊗ₜ[ℕ] m) = Pf.mk m a
  rw [pfTToPf_tmul_of 1 a (pnatInv a) (pnatInv_mul a) m, one_smul]

/-- ★★**往復 2** —— `ℚ≥0 ⊗ M → M^pf → ℚ≥0 ⊗ M` は恒等。 -/
theorem pfToPfT_pfTToPf (x : PfT M) : pfToPfT (pfTToPf x) = x := by
  induction x using TensorProduct.induction_on with
  | zero => rw [map_zero, map_zero]
  | tmul q m =>
      rw [pfTToPf_tmul]
      show pnatInv (nnratDen q) ⊗ₜ[ℕ] (q.num • m) = q ⊗ₜ[ℕ] m
      rw [← TensorProduct.smul_tmul]
      congr 1
      have hq : ((q.den : ℕ) : ℚ≥0) ≠ 0 := by exact_mod_cast q.den_pos.ne'
      rw [nsmul_eq_mul, pnatInv]
      show ((q.num : ℕ) : ℚ≥0) * (((q.den : ℕ)) : ℚ≥0)⁻¹ = q
      rw [nnrat_num_eq q]
      field_simp
  | add x y hx hy => rw [map_add, map_add, hx, hy]

/-- ★★★**`M^pf ≃ ℚ≥0 ⊗_ℕ M`** —— 商模型とテンソル模型の一致。 -/
noncomputable def pfEquivPfT : Pf M ≃+ PfT M where
  toFun := pfToPfT
  invFun := pfTToPf
  left_inv := pfTToPf_pfToPfT
  right_inv := pfToPfT_pfTToPf
  map_add' x y := map_add pfToPfT x y

/-- ★同一視は**自然**である。 -/
theorem pfToPfT_natural (f : M →+ N) (x : Pf M) :
    pfToPfT (Pf.map f x) = scMap (S := ℚ≥0) f (pfToPfT x) := by
  induction x using Pf.inductionOn with | _ m a =>
  show pfToPfT (Pf.mk (f m) a) = scMap (S := ℚ≥0) f (pnatInv a ⊗ₜ[ℕ] m)
  rw [scMap_tmul]
  rfl

/-! ## ★6. 帰結 —— 条件 (a) が `S = ℚ≥0` では**証明できる** -/

/-- ★★★★★**`ℚ≥0 ⊗_ℕ −` は単射性を保つ**。

★★平坦性を直接示したのではなく、商模型 `M^pf` へ移して
`Pf.map_injective`(5 行)に帰着させたものである。 -/
theorem scMap_injective_nnrat {f : M →+ N} (hf : Function.Injective f) :
    Function.Injective (scMap (S := ℚ≥0) f) := by
  intro x y h
  have hx := pfToPfT_pfTToPf x
  have hy := pfToPfT_pfTToPf y
  rw [← hx, ← hy, ← pfToPfT_natural, ← pfToPfT_natural] at h
  have h2 : Pf.map f (pfTToPf x) = Pf.map f (pfTToPf y) :=
    pfEquivPfT.injective h
  rw [← hx, ← hy, Pf.map_injective hf h2]

theorem isSharp_of_addEquiv {X Y : Type*} [AddCommMonoid X] [AddCommMonoid Y]
    (e : X ≃+ Y) (h : IsSharp X) : IsSharp Y := by
  intro a ha
  have h1 : IsAddUnit (e.symm a) := ha.map e.symm.toAddMonoidHom
  have h2 : e.symm a = 0 := h _ h1
  have := congrArg e h2
  simpa using this

/-- ★★**`ℚ≥0 ⊗_ℕ M` は `M` が sharp なら sharp**。 -/
theorem isSharp_pfT (h : IsSharp M) : IsSharp (PfT M) :=
  isSharp_of_addEquiv pfEquivPfT (Pf.isSharp_pf h)

/-- ★★★★★★**`ℚ≥0 ⊗_ℕ −` は characteristically injective を保つ** ——
`Def24Rlf.lean` の `phiScOn` が仮定に置いていた `hcharInj` の、`S = ℚ≥0` の場合。 -/
theorem isCharacteristicallyInjective_scMap_nnrat (hN : IsSharp N)
    {f : M →+ N} (hf : Function.Injective f) :
    IsCharacteristicallyInjective (scMap (S := ℚ≥0) f) :=
  ⟨scMap_injective_nnrat hf,
    charMap_injective_of_sharp (isSharp_pfT hN) (scMap_injective_nnrat hf)⟩

/-! ## ★7. `Φ^pf : MonoidOn 𝒟` を**仮定なしで** -/

theorem isCancelAdd_of_addEquiv {X Y : Type*} [AddCommMonoid X] [AddCommMonoid Y]
    (e : X ≃+ Y) [IsCancelAdd X] : IsCancelAdd Y := by
  haveI : IsLeftCancelAdd Y := by
    refine ⟨fun a b c h => ?_⟩
    have h0 : a + b = a + c := h
    have h2 : e.symm a + e.symm b = e.symm a + e.symm c := by
      rw [← map_add, ← map_add, h0]
    have := add_left_cancel h2
    simpa using congrArg e this
  haveI : IsRightCancelAdd Y := by
    refine ⟨fun a b c h => ?_⟩
    have h0 : b + a = c + a := h
    have h2 : e.symm b + e.symm a = e.symm c + e.symm a := by
      rw [← map_add, ← map_add, h0]
    have := add_right_cancel h2
    simpa using congrArg e this
  exact ⟨⟩

/-- ★★**`ℚ≥0 ⊗_ℕ M` は `M` が integral なら integral**。

★`Prop53Rlf.lean` の `scModel` が仮定に置いていた `hint` の、`S = ℚ≥0` の場合。 -/
theorem isIntegralMonoid_pfT (h : IsIntegralMonoid M) : IsIntegralMonoid (PfT M) := by
  letI := isCancelAdd_of_isIntegralMonoid M h
  letI : IsCancelAdd (Pf M) := Pf.isCancelAdd'
  letI : IsCancelAdd (PfT M) := isCancelAdd_of_addEquiv (pfEquivPfT (M := M))
  exact isIntegralMonoid_of_isCancelAdd (PfT M)

variable {D : Type u} [Category.{v} D]

/-- ★★★★★★★**`Φ^pf = ℚ≥0 ⊗_ℕ Φ` は仮定なしで `MonoidOn 𝒟` になる**。

★`Def24Rlf.lean` の `phiScOn` は条件 (a) を仮定 `hcharInj` として受け取っていたが、
`S = ℚ≥0` ではそれが**証明できる**。
★残る `rlf-flat` は **`S = ℝ≥0` の場合だけ**である。 -/
noncomputable def phiPfTOn (Φ : MonoidOn.{v, u, w} D)
    (hsharp : ∀ A : D, IsSharp (Φ.val A)) : MonoidOn.{v, u, w} D :=
  phiScOn NNRat Φ (fun {_ B} α =>
    isCharacteristicallyInjective_scMap_nnrat (hsharp B) (Φ.map_injective α))

/-- ★divisorial なら各 `Φ(A)` は sharp なので、仮定は要らない。 -/
noncomputable def phiPfTOnOfDivisorial (Φ : MonoidOn.{v, u, w} D)
    (hdiv : Φ.IsDivisorialOn) : MonoidOn.{v, u, w} D :=
  phiPfTOn Φ (fun A => (hdiv A).2)

/-! ### ★出典の紐付け -/

/-- ★locator —— `M^pf` の 2 つの模型の一致。 -/
def pfEquivPfT.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 19, item := "Definition 1.1, (ii) — Φ^pf(商模型とテンソル模型の一致)",
    sectionId := "frdi-def-1-1-ii" }

/-- ★locator —— 条件 (a) の `S = ℚ≥0` の場合。 -/
def phiPfTOn.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 19, item := "Definition 1.1, (ii) — Φ^pf は monoid on 𝒟",
    sectionId := "frdi-def-1-1-ii" }

end ABC3.Found.FrdI
