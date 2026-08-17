import ABC3.Found.Arakelov.PicUnitFree

/-!
# Arakelov (B1) 第 32 ブロック —— **押し出しの `μ` は切断の上で恒等である**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★`δ` の定義式の残り 1 つ

    δ P Q = homEquiv.symm ((unit P ⊗ₘ unit Q) ≫ μ)

★第 31 ブロックで `unit` を生成元の上で書き下した。
★★本ブロックで `μ` を**切断の上で**書き下す。

## ★★★★μ は「何もしない」

押し出しの切断は `(f_* M)(V) = M(f⁻¹V)` であり、
テンソルは切断ごとなので、`μ` は

    m ⊗ₜ n ↦ m ⊗ₜ n

である。★★ただし**係数環が変わる**(`R_Y(V) → R_X(f⁻¹V)`)ので、
`rfl` では通らない——`ModuleCat.restrictScalars_μ_tmul` を 1 回当てる。

★★★これは第 18 ブロック(`resScalars` の lax 構造)の recipe そのものである。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `pushforward_μ_app_tmul` | ★★★★**`μ` は純テンソルを純テンソルに送る** |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace Limits

variable {X Y : Scheme.{u}} (f : X ⟶ Y)

/-- ★押し出しの lax monoidal 構造(第 19 ブロック)をインスタンスとして置く。 -/
noncomputable instance pushforwardLax :
    (PresheafOfModules.pushforward.{u} (pullbackPhi f)).LaxMonoidal :=
  pushCRLax (Opens.map f.base) X.presheaf Y.presheaf f.c

/-- ★★★★**押し出しの `μ` は切断の上で純テンソルを純テンソルに送る**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★これで `δ` の定義式の右辺(`(unit ⊗ₘ unit) ≫ μ`)が
**切断の上で計算できる**ようになった。

★係数環が `R_Y(V)` から `R_X(f⁻¹V)` へ変わるので `rfl` では通らない。
`ModuleCat.restrictScalars_μ_tmul` を 1 回当てる——第 18 ブロックの recipe である。 -/
theorem pushforward_μ_app_tmul (M N : X.PresheafOfModules) (V : Y.Opens)
    (m : ((PresheafOfModules.pushforward (pullbackPhi f)).obj M).obj (op V))
    (n : ((PresheafOfModules.pushforward (pullbackPhi f)).obj N).obj (op V)) :
    (Functor.LaxMonoidal.μ (PresheafOfModules.pushforward (pullbackPhi f)) M N).app (op V)
        (m ⊗ₜ n)
      = (show ((PresheafOfModules.pushforward (pullbackPhi f)).obj (M ⊗ N)).obj (op V) from
          m ⊗ₜ n) := by
  erw [ModuleCat.restrictScalars_μ_tmul]
  rfl

/-! ## ★出典の紐付け(`.src`) -/

def pushforward_μ_app_tmul.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——押し出しの μ が切断の上で恒等であること)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
