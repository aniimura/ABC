import ABC3.Found.Arakelov.UltraPoint

/-!
# Arakelov (C2) 第 83 ブロック —— **★★★chart 間の翻訳は係数環に自然**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★超積で困るのはここだった

超積体の点は chart `U` から作った。★しかし**極限点は `U` の外に出る**
——`U` の閉包に入るだけである。★★したがって極限の近傍を測る chart `V` は `U` と違う。

★★★**`U` の座標と `V` の座標を繋ぐ必要がある。**

## ★★★★★★鍵は「係数環に自然」であること

`φ : Γ(X,U) → R` が定める点 `Spec R ⟶ X` が `V` に入るなら、
`V` 座標 `Ψ(φ) : Γ(X,V) → R` が定まる。★このとき

    Ψ(ρ ∘ φ) = ρ ∘ Ψ(φ)      (ρ : R → R' 環準同型)

★★これは**開埋め込みへの持ち上げの一意性**だけから出る。

★★★これで「超積の座標 = 座標の超積」が言える——
`ρ` として **芽を取る準同型** `(Arc X → ℂ) → *ℂ` を入れればよい。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `liftV` / `liftV_fac` | ★開埋め込みへの持ち上げ |
| `liftV_comp` | ★★持ち上げは合成と両立(一意性) |
| `chartPoint` | ★chart の環準同型が定める点 |
| `chartHom` | ★★`U` 座標 → `V` 座標 |
| `chartHom_natural` | ★★★**係数環への自然性** |
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory

variable {X : Scheme.{0}}

/-! ## ★開埋め込みへの持ち上げ -/

/-- ★開埋め込みへの持ち上げ。 -/
noncomputable def liftV {W : Scheme.{0}} (V : X.Opens) (g : W ⟶ X)
    (h : Set.range g.base ⊆ (V : Set X)) : W ⟶ V.toScheme :=
  IsOpenImmersion.lift V.ι g (by rw [Scheme.Opens.range_ι]; exact h)

theorem liftV_fac {W : Scheme.{0}} (V : X.Opens) (g : W ⟶ X)
    (h : Set.range g.base ⊆ (V : Set X)) : liftV V g h ≫ V.ι = g :=
  IsOpenImmersion.lift_fac _ _ _

/-- ★★持ち上げは合成と両立する(一意性から)。 -/
theorem liftV_comp {W W' : Scheme.{0}} (V : X.Opens) (g : W ⟶ X)
    (h : Set.range g.base ⊆ (V : Set X)) (e : W' ⟶ W)
    (h' : Set.range (e ≫ g).base ⊆ (V : Set X)) :
    liftV V (e ≫ g) h' = e ≫ liftV V g h := by
  apply (cancel_mono V.ι).1
  rw [liftV_fac, Category.assoc, liftV_fac]

/-! ## ★★chart 間の翻訳 -/

/-- ★`U` chart の環準同型が定める `Spec R ⟶ X`。 -/
noncomputable def chartPoint {R : CommRingCat.{0}} (U : X.affineOpens) (φ : Γ(X, U.1) ⟶ R) :
    Spec R ⟶ X :=
  Spec.map φ ≫ U.2.isoSpec.inv ≫ U.1.ι

theorem chartPoint_comp {R R' : CommRingCat.{0}} (U : X.affineOpens) (φ : Γ(X, U.1) ⟶ R)
    (ρ : R ⟶ R') : chartPoint U (φ ≫ ρ) = Spec.map ρ ≫ chartPoint U φ := by
  rw [chartPoint, chartPoint, Spec.map_comp, Category.assoc]

/-- ★★`U` chart の環準同型を `V` chart の環準同型に翻訳する。 -/
noncomputable def chartHom {R : CommRingCat.{0}} (U V : X.affineOpens) (φ : Γ(X, U.1) ⟶ R)
    (h : Set.range (chartPoint U φ).base ⊆ (V.1 : Set X)) : Γ(X, V.1) ⟶ R :=
  Spec.preimage (liftV V.1 (chartPoint U φ) h ≫ V.2.isoSpec.hom)

/-- ★★★**翻訳は係数環に自然である**——これが超積の要。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★「超積の座標 = 座標の超積」がこれから出る。 -/
theorem chartHom_natural {R R' : CommRingCat.{0}} (U V : X.affineOpens)
    (φ : Γ(X, U.1) ⟶ R) (ρ : R ⟶ R')
    (h : Set.range (chartPoint U φ).base ⊆ (V.1 : Set X))
    (h' : Set.range (chartPoint U (φ ≫ ρ)).base ⊆ (V.1 : Set X)) :
    chartHom U V (φ ≫ ρ) h' = chartHom U V φ h ≫ ρ := by
  have hcp : chartPoint U (φ ≫ ρ) = Spec.map ρ ≫ chartPoint U φ := chartPoint_comp U φ ρ
  unfold chartHom
  have hlift : liftV V.1 (chartPoint U (φ ≫ ρ)) h'
      = Spec.map ρ ≫ liftV V.1 (chartPoint U φ) h := by
    apply (cancel_mono V.1.ι).1
    rw [liftV_fac, Category.assoc, liftV_fac, hcp]
  rw [hlift, Category.assoc, Spec.preimage_comp, Spec.preimage_map]

/-! ## ★出典の紐付け(`.src`) -/

def chartHom_natural.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 C——chart 間の座標変換の係数環自然性)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
