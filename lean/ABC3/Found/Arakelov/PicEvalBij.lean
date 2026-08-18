import ABC3.Found.Arakelov.PicDualModules

/-!
# Arakelov (B2) 第 181 ブロック —— **自明な所では評価射は全単射**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★評価射が同型になる理由

`F|_V ≅ 𝟙_` なる `V` の上では

    F(V) ⊗ F^∨(V) → 𝒪(V)、  x ⊗ φ ↦ φ(x)

は全単射である。★★理由は**双対の切断がすべて `c • e.hom` の形**だから
(第 171 の `bijective_dual_smul`)——したがって

    x ⊗ φ = x ⊗ (c • e.hom) = (c • x) ⊗ e.hom

と**生成元 1 個**に潰れ、逆写像は `a ↦ e⁻¹(a) ⊗ e.hom` で与えられる。

## ★★逆写像を手で書くのが一番短い

同型の合成(`F(V) ⊗ 𝒪(V) ≅ F(V) ≅ 𝒪(V)`)で組むと**係数環の綴りが 3 通り**
現れる([[ring-instance-two-paths]])。★逆写像を直に書き、
`TensorProduct.induction_on` で 3 場合を潰す方が短い。

## ★★本ブロックで取れるもの

| 定義・定理 | 内容 |
|---|---|
| `eSec` / `gVal` | ★型の橋 |
| `dual_smul_surj` | ★双対の切断はすべて `c • e.hom` |
| `evHom_tmul_hom` / `evHom_tmul_smul` | ★評価射の値(`rfl` と `map_smul`) |
| `evInvFun` | ★★逆写像 |
| `bijective_evHom_app` | ★★★★★**自明な所では評価射は全単射** |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace

variable {X : Scheme.{u}} (F : X.PresheafOfModules) (V : X.Opens)
  (e : (restrictPresheafFunctor X V).obj F ≅ 𝟙_ (PresheafModulesOn X V))

/-- ★**`e.hom` を双対の切断として見る**。 -/
abbrev eSec : ((dualPresheaf F).obj (op V) : Type u) := e.hom

/-- ★**逆向きの橋**(`fVal` の逆)。 -/
def gVal (W : (X.Opens)ᵒᵖ)
    (x : (((restrictPresheafFunctor X W.unop).obj F).obj (op (Over.mk (𝟙 W.unop))) : Type u)) :
    (F.obj W : Type u) := x

/-- ★**双対の切断はすべて `c • e.hom` の形である**。 -/
theorem dual_smul_surj (φ : ((dualPresheaf F).obj (op V) : Type u)) :
    ∃ c : (Γ(X, V) : Type u), c • eSec F V e = φ :=
  (bijective_dual_smul F V e (Over.mk (𝟙 V))).2 φ

/-- ★**評価射の値(生成元)**(`rfl`)。 -/
theorem evHom_tmul_hom (x : (F.obj (op V) : Type u)) :
    ((evHom F).app (op V)).hom
        (TensorProduct.tmul (Γ(X, V) : Type u) x (eSec F V e))
      = ((e.hom.app (op (Over.mk (𝟙 V)))).hom (fVal F (op V) x) : (Γ(X, V) : Type u)) := rfl

/-- ★**評価射の値(係数つき)**。 -/
theorem evHom_tmul_smul (x : (F.obj (op V) : Type u)) (c : (Γ(X, V) : Type u)) :
    ((evHom F).app (op V)).hom
        (TensorProduct.tmul (Γ(X, V) : Type u) x (c • eSec F V e))
      = rVal (X := X) (op V) c • ((e.hom.app (op (Over.mk (𝟙 V)))).hom (fVal F (op V) x)) :=
  (evBil F (op V) x).map_smul c (eSec F V e)

/-- ★同型の 2 法則(値の側)。 -/
theorem eA_eB (x : (F.obj (op V) : Type u)) :
    (e.inv.app (op (Over.mk (𝟙 V)))).hom
        ((e.hom.app (op (Over.mk (𝟙 V)))).hom (fVal F (op V) x)) = fVal F (op V) x :=
  congrArg (fun (m : _ ⟶ _) => ((m.app (op (Over.mk (𝟙 V)))).hom (fVal F (op V) x))) e.hom_inv_id

theorem eB_eA (a : (Γ(X, V) : Type u)) :
    (e.hom.app (op (Over.mk (𝟙 V)))).hom
        ((e.inv.app (op (Over.mk (𝟙 V)))).hom a) = a :=
  congrArg (fun (m : _ ⟶ _) => ((m.app (op (Over.mk (𝟙 V)))).hom a)) e.inv_hom_id

/-- ★★**評価射の逆写像**——`a ↦ e⁻¹(a) ⊗ e.hom`。 -/
noncomputable def evInvFun (a : (Γ(X, V) : Type u)) :
    (((F ⊗ dualPresheaf F).obj (op V)) : Type u) :=
  TensorProduct.tmul (Γ(X, V) : Type u)
    (gVal F (op V) ((e.inv.app (op (Over.mk (𝟙 V)))).hom a)) (eSec F V e)

theorem evInvFun_add (a b : (Γ(X, V) : Type u)) :
    evInvFun F V e (a + b) = evInvFun F V e a + evInvFun F V e b :=
  (congrArg (fun y => TensorProduct.tmul (Γ(X, V) : Type u) (gVal F (op V) y) (eSec F V e))
      (map_add (e.inv.app (op (Over.mk (𝟙 V)))).hom a b)).trans
    (TensorProduct.add_tmul _ _ _)

/-- ★評価射は逆写像の左逆である。 -/
theorem evHom_evInvFun (a : (Γ(X, V) : Type u)) :
    ((evHom F).app (op V)).hom (evInvFun F V e a) = a :=
  (evHom_tmul_hom F V e _).trans (eB_eA F V e a)

/-- ★★逆写像は評価射の左逆である——`induction_on` で 3 場合。 -/
theorem evInvFun_evHom (t : (((F ⊗ dualPresheaf F).obj (op V)) : Type u)) :
    evInvFun F V e (((evHom F).app (op V)).hom t) = t := by
  induction t using TensorProduct.induction_on with
  | zero =>
      refine Eq.trans (congrArg (evInvFun F V e) (map_zero ((evHom F).app (op V)).hom)) ?_
      refine Eq.trans (congrArg
        (fun y => TensorProduct.tmul (Γ(X, V) : Type u) (gVal F (op V) y) (eSec F V e))
        (map_zero (e.inv.app (op (Over.mk (𝟙 V)))).hom)) ?_
      exact TensorProduct.zero_tmul _ _
  | tmul x φ =>
      obtain ⟨c, rfl⟩ := dual_smul_surj F V e φ
      refine Eq.trans (congrArg (evInvFun F V e) (evHom_tmul_smul F V e x c)) ?_
      refine Eq.trans (congrArg
        (fun y => TensorProduct.tmul (Γ(X, V) : Type u) (gVal F (op V) y) (eSec F V e))
        ((e.inv.app (op (Over.mk (𝟙 V)))).hom.map_smul (rVal (X := X) (op V) c) _)) ?_
      refine Eq.trans (congrArg
        (fun y => TensorProduct.tmul (Γ(X, V) : Type u)
          (gVal F (op V) (rVal (X := X) (op V) c • y)) (eSec F V e))
        (eA_eB F V e x)) ?_
      exact (TensorProduct.smul_tmul' c x (eSec F V e)).trans
        (TensorProduct.tmul_smul c x (eSec F V e)).symm
  | add a b ha hb =>
      refine Eq.trans (congrArg (evInvFun F V e) (map_add ((evHom F).app (op V)).hom a b)) ?_
      exact (evInvFun_add F V e _ _).trans (congrArg₂ (· + ·) ha hb)
include e in

/-- ★★★★★**自明な所では評価射は全単射である**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★これで「評価射は局所全単射」が言え、層化で同型になる。 -/
theorem bijective_evHom_app : Function.Bijective ((evHom F).app (op V)).hom :=
  Function.bijective_iff_has_inverse.2
    ⟨evInvFun F V e, evInvFun_evHom F V e, evHom_evInvFun F V e⟩

/-! ## ★出典の紐付け(`.src`) -/

def bijective_evHom_app.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——自明な所では評価射は全単射)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
