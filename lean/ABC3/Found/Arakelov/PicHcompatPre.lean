import ABC3.Found.Arakelov.PicSpecMapLE

/-!
# Arakelov (B2) 第 229 ブロック —— ★★★★★★**`hcompat` が完成した**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★7 つの摩擦を越えて到達した

    hB.fromSpec.app B ≫ 転送 ≫ ΓSpecIso.hom ≫ f.appLE B A i ≫ ΓSpecIso.inv ≫ 制限
      = f.app B ≫ hA.fromSpec.appLE _ _ e₂

★これが第 221 が受けた仮定 `hcompat` の中身である。

## ★★機構は 2 本の合流

| 段 | ブロック |
|---|---|
| 可換正方形(座標を揃えた形) | ★第 227 |
| `Spec.map` の `appLE` を `ΓSpecIso` で書く | ★第 228 |

★★`rw [← specMap_appLE]` で右半分を `Spec.map` の形に戻し、第 227 を当てるだけ。

## ★★★この 1 欄で越えた摩擦の一覧

| # | 症状 | 逃げ道 | ブロック |
|---|---|---|---|
| 1 | `whnf` timeout | 抽象の側で書く | 219 |
| 2 | 始域が違う | 分解を挟む | 215 |
| 3 | `simp` の過剰正規化 | 元のレベルで書く | 217 |
| 3' | ★`rw` が結合で詰まる | ★**`simp` に任せる** | ★228 |
| 4 | 依存位置 | grep し直す / `congr_app` | 222 / 225 |
| 5 | 転送の伝播 | 同型であることを使う | 223 |
| 6 | `rw` が証明部分まで照合 | 証明を引数に出す | 224 |
| 7 | 座標系が合わない | `appLE` で依存しない座標 | 226–227 |

★★★★★★**8 通りの逃げ道**が要った。これが `ofDivisor_pullback` 1 欄の実態である。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `hcompat_pre` | ★★★★★★**`hcompat`** |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite TopologicalSpace

variable {X Y : Scheme.{u}} (f : X ⟶ Y)

/-- ★★★★★★**`hcompat`**(座標を `fromSpec⁻¹` に取った形)。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★第 227(可換正方形)と第 228(`Spec.map` の `appLE`)の合流である。 -/
theorem hcompat_pre {A : X.Opens} {B : Y.Opens}
    (hA : IsAffineOpen A) (hB : IsAffineOpen B) (i : A ≤ f ⁻¹ᵁ B)
    (e1 : (hA.fromSpec ⁻¹ᵁ A) ≤ (Spec.map (f.appLE B A i) ≫ hB.fromSpec) ⁻¹ᵁ B)
    (e2 : (hA.fromSpec ⁻¹ᵁ A) ≤ (hA.fromSpec ≫ f) ⁻¹ᵁ B)
    (e3 : (hA.fromSpec ⁻¹ᵁ A) ≤ (Spec.map (f.appLE B A i)) ⁻¹ᵁ ⊤) :
    hB.fromSpec.app B
        ≫ (Spec Γ(Y, B)).presheaf.map (eqToHom hB.fromSpec_preimage_self.symm).op
        ≫ (Scheme.ΓSpecIso (Γ(Y, B))).hom ≫ f.appLE B A i
        ≫ (Scheme.ΓSpecIso (Γ(X, A))).inv
        ≫ (Spec Γ(X, A)).presheaf.map (homOfLE le_top).op
      = f.app B ≫ hA.fromSpec.appLE _ _ e2 := by
  rw [← specMap_appLE (f.appLE B A i) _ e3]
  exact hcompat_core f hA hB i e1 e2


/-! ## ★出典の紐付け(`.src`) -/

def hcompat_pre.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——hcompat)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
