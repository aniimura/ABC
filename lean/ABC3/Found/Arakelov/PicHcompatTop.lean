import ABC3.Found.Arakelov.PicPreCoord

/-!
# Arakelov (B2) 第 231 ブロック —— ★★★★★★**`hcompat` を `⊤` 座標で**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★§9-263 の教訓が**自分の補題にも当たった**

§9-262 では「第 229 は `fromSpec⁻¹` 座標だから、第 211/220 も揃える」と考えた。
★★しかし**第 226(`square_appLE`)は座標 `W` をパラメータに取ってある**——
★★★**`W := ⊤` と取り直すだけ**で、第 210/211/220 と座標が揃う。

★§9-263 で「mathlib の補題が座標をパラメータに取っている」と気づいたが、
**自前の第 226 も同じだった**——**自分が第 227 で `fromSpec⁻¹A` に固定していた**。

## ★★これで座標変更のコストが 0 になった

| 補題 | 座標 | 変更コスト |
|---|---|---|
| 第 226 `square_appLE` | ★パラメータ | ★**0** |
| 第 228 `specMap_appLE` | ★パラメータ(`W`) | ★**0** |
| 第 219 `ideal_of_comap_eq` | ★抽象 | ★**0** |
| 第 218 `surjective_pullIdealHom` | ★パラメータ | ★**0** |

★★★★**パラメータに取っておいた 4 本すべてが、座標変更を無料にした**。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `hcompat_top` | ★★★★★`⊤` 座標の可換正方形 |
| `hcompat_top_gamma` | ★★★★★★**`⊤` 座標の `hcompat`** |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite TopologicalSpace

variable {X Y : Scheme.{u}} (f : X ⟶ Y)

theorem hcompat_top {A : X.Opens} {B : Y.Opens}
    (hA : IsAffineOpen A) (hB : IsAffineOpen B) (i : A ≤ f ⁻¹ᵁ B)
    (e1 : (⊤ : (Spec Γ(X, A)).Opens) ≤ (Spec.map (f.appLE B A i) ≫ hB.fromSpec) ⁻¹ᵁ B)
    (e2 : (⊤ : (Spec Γ(X, A)).Opens) ≤ (hA.fromSpec ≫ f) ⁻¹ᵁ B) :
    hB.fromSpec.app B ≫ (Spec.map (f.appLE B A i)).appLE _ _ e1
      = f.app B ≫ hA.fromSpec.appLE _ _ e2 := by
  rw [← Scheme.Hom.comp_appLE, ← Scheme.Hom.comp_appLE]
  exact square_appLE f hA hB i ⊤ e1 e2




/-- ★★★★★★`hcompat`(`⊤` 座標、`ΓSpecIso` で書いた形)。 -/
theorem hcompat_top_gamma {A : X.Opens} {B : Y.Opens}
    (hA : IsAffineOpen A) (hB : IsAffineOpen B) (i : A ≤ f ⁻¹ᵁ B)
    (e1 : (⊤ : (Spec Γ(X, A)).Opens) ≤ (Spec.map (f.appLE B A i) ≫ hB.fromSpec) ⁻¹ᵁ B)
    (e2 : (⊤ : (Spec Γ(X, A)).Opens) ≤ (hA.fromSpec ≫ f) ⁻¹ᵁ B)
    (e3 : (⊤ : (Spec Γ(X, A)).Opens) ≤ (Spec.map (f.appLE B A i)) ⁻¹ᵁ ⊤) :
    hB.fromSpec.app B
        ≫ (Spec Γ(Y, B)).presheaf.map (eqToHom hB.fromSpec_preimage_self.symm).op
        ≫ (Scheme.ΓSpecIso (Γ(Y, B))).hom ≫ f.appLE B A i
        ≫ (Scheme.ΓSpecIso (Γ(X, A))).inv
        ≫ (Spec Γ(X, A)).presheaf.map (homOfLE le_top).op
      = f.app B ≫ hA.fromSpec.appLE _ _ e2 := by
  rw [← specMap_appLE (f.appLE B A i) _ e3]
  exact hcompat_top f hA hB i e1 e2


/-! ## ★出典の紐付け(`.src`) -/

def hcompat_top_gamma.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——hcompat を ⊤ 座標で)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
