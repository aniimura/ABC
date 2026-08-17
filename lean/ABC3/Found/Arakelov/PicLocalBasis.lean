import ABC3.Found.Arakelov.PicLocality

/-!
# Arakelov (B1) 第 9 ブロック —— **局所基底からの局所単射性**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★`Over` を経由しない

当初は「制限(`Over V`)+ 制限のモノイダル性 + 篩の往復」で
局所自明化を接続するつもりだった。★★**それは要らなかった。**

★★★実際に要るのは **`V` の上での基底 1 本だけ**である:

    x = x' ⊗ s,  y = y' ⊗ s   (s は M(V) の基底)
    ⟹ f(x') ⊗ s = f(y') ⊗ s ⟹ f(x') = f(y')
    ⟹ equalizerSieve x' y' ∈ J V     (f の局所単射性)
    ⟹ equalizerSieve x y ⊇ それ       (`tensorObj_map_tmul` だけ)

★`W ⟶ V` の上で基底であることは**要らない**——制限は `tensorObj_map_tmul` で
そのまま計算できるからである。

## ★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `tmul_basis_eq` | 階数 1 自由なら `z = φ z ⊗ₜ e 1` |
| `isLocallyInjective_whiskerRight_of_basis` | ★**局所基底があれば局所単射性はテンソルで保たれる** |
-/

universe u

namespace ABC3.Found.Arakelov

open CategoryTheory MonoidalCategory Opposite TensorProduct

variable {C : Type u} [Category.{u} C] (J : GrothendieckTopology C)
  {R : Cᵒᵖ ⥤ CommRingCat.{u}}

/-! ## ★★階数 1 自由加群とのテンソル積 -/

section Module

variable {A : Type u} [CommRing A] {N L : Type u} [AddCommGroup N] [Module A N]
  [AddCommGroup L] [Module A L]

/-- ★`L ≃ A` のとき、`N ⊗ L ≃ N`。 -/
noncomputable def tensorRankOne (e : A ≃ₗ[A] L) : N ⊗[A] L ≃ₗ[A] N :=
  (TensorProduct.congr (LinearEquiv.refl A N) e.symm).trans (TensorProduct.rid A N)

/-- ★★**階数 1 自由なら、すべての元は `x ⊗ e 1` の形である**。 -/
theorem tmul_basis_eq (e : A ≃ₗ[A] L) (z : N ⊗[A] L) :
    tensorRankOne (N := N) e z ⊗ₜ[A] e 1 = z := by
  induction z using TensorProduct.induction_on with
  | zero => simp [tensorRankOne]
  | tmul n l =>
      show (TensorProduct.rid A N) ((TensorProduct.congr (LinearEquiv.refl A N) e.symm)
        (n ⊗ₜ l)) ⊗ₜ[A] e 1 = n ⊗ₜ l
      simp only [TensorProduct.congr_tmul, LinearEquiv.refl_apply, TensorProduct.rid_tmul]
      rw [TensorProduct.smul_tmul, ← map_smul, smul_eq_mul, mul_one,
        LinearEquiv.apply_symm_apply]
  | add z w hz hw =>
      rw [map_add, TensorProduct.add_tmul, hz, hw]

/-- ★★**基底とのテンソルは単射的**——`x ⊗ e 1 = y ⊗ e 1 → x = y`。 -/
theorem tmul_basis_injective (e : A ≃ₗ[A] L) {x y : N}
    (h : x ⊗ₜ[A] e 1 = y ⊗ₜ[A] e 1) : x = y := by
  have := congrArg (tensorRankOne (N := N) e) h
  simpa [tensorRankOne] using this

end Module

/-! ## ★★★★局所基底からの局所単射性 -/

/-- ★★★★★**局所基底があれば、局所単射性はテンソルで保たれる**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★これが結合律の残り 1 条である。仮定は
「各 `U` に被覆篩があり、その各員 `V` の上で `M(V)` が階数 1 自由」——
すなわち **`M` が可逆層である**ことの言い換えである。 -/
theorem isLocallyInjective_whiskerRight_of_basis
    {P Q : PresheafOfModules.{u} (R ⋙ forget₂ CommRingCat.{u} RingCat.{u})}
    (f : P ⟶ Q) [PresheafOfModules.IsLocallyInjective J f]
    (M : PresheafOfModules.{u} (R ⋙ forget₂ CommRingCat.{u} RingCat.{u}))
    (htriv : ∀ U : C, ∃ S : Sieve U, S ∈ J U ∧ ∀ ⦃V : C⦄ (i : V ⟶ U), S i →
      Nonempty ((R.obj (op V) : Type u) ≃ₗ[(R.obj (op V) : Type u)] M.obj (op V))) :
    PresheafOfModules.IsLocallyInjective J (f ▷ M) := by
  refine isLocallyInjective_of_cover J _ ?_
  intro U x y hxy
  obtain ⟨S, hS, htv⟩ := htriv U
  refine ⟨S, hS, ?_⟩
  intro V i hi
  obtain ⟨e⟩ := htv i hi
  -- ★制限した切断を基底で分解する
  set xV := (P ⊗ M).map i.op x with hxVdef
  set yV := (P ⊗ M).map i.op y with hyVdef
  set x' := tensorRankOne (N := (P.obj (op V) : Type u)) e xV with hx'
  set y' := tensorRankOne (N := (P.obj (op V) : Type u)) e yV with hy'
  have hxV : x' ⊗ₜ e 1 = xV := tmul_basis_eq e xV
  have hyV : y' ⊗ₜ e 1 = yV := tmul_basis_eq e yV
  -- ★制限しても像は等しい(自然性)
  have hnat : (f ▷ M).app (op V) xV = (f ▷ M).app (op V) yV := by
    rw [hxVdef, hyVdef]
    show ((PresheafOfModules.toPresheaf _).map (f ▷ M)).app (op V)
        (((PresheafOfModules.toPresheaf _).obj (P ⊗ M)).map i.op x)
      = ((PresheafOfModules.toPresheaf _).map (f ▷ M)).app (op V)
        (((PresheafOfModules.toPresheaf _).obj (P ⊗ M)).map i.op y)
    rw [NatTrans.naturality_apply, NatTrans.naturality_apply, hxy]
  -- ★基底を落として `f` の等式にする
  have hf : f.app (op V) x' = f.app (op V) y' := by
    refine tmul_basis_injective e ?_
    have key : (f ▷ M).app (op V) (x' ⊗ₜ e 1) = (f ▷ M).app (op V) (y' ⊗ₜ e 1) := by
      rw [hxV, hyV]
      exact hnat
    exact key
  -- ★`f` の局所単射性を使う
  refine J.superset_covering ?_
    (Presheaf.equalizerSieve_mem J ((PresheafOfModules.toPresheaf _).map f) x' y' hf)
  intro W j hj
  show (P ⊗ M).map j.op xV = (P ⊗ M).map j.op yV
  rw [← hxV, ← hyV]
  have ht1 := PresheafOfModules.Monoidal.tensorObj_map_tmul (M₁ := P) (M₂ := M) j.op x' (e 1)
  have ht2 := PresheafOfModules.Monoidal.tensorObj_map_tmul (M₁ := P) (M₂ := M) j.op y' (e 1)
  have hjj : P.map j.op x' = P.map j.op y' := hj
  refine ht1.trans (Eq.trans ?_ ht2.symm)
  rw [hjj]

/-! ## ★出典の紐付け(`.src`) -/

def tmul_basis_eq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——階数 1 自由加群とのテンソル積の分解)",
    sectionId := "genell-def-1-1-i" }

def isLocallyInjective_whiskerRight_of_basis.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——局所基底からの局所単射性)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
