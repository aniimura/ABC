import ABC3.Found.GaloisRep.CurveHomVal
import ABC3.Found.GaloisRep.TranslateAutAll

/-!
# Galois (G5) 第 168 ブロック —— **★★★★★★★`μ f_P` は平行移動で不変**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★`[n] ∘ τ_T = [n]` を Lean の言葉にした

`T ∈ E[n]` に対する平行移動 `τ_T`(第 124 で関数体の `F`-自己同型として構成済)は

    τ_T ∘ μ = μ         (μ = [n]^*)

を満たす。★これが `e_v` の一定性の**唯一の源**である。

### ★★★★★機構は 4 段

1. **`τ_T` は生成点を `生成点 + T` に送る**(`autFF_generic`)——
   mathlib の `Point.map`(**加法群の準同型**)に第 124 の
   `τ(coordX) = translateX`、`τ(coordY) = translateY` を当てる。
   ★`translateX`・`translateY` が「生成点 + 定数点」であることは
   `Point.add_some` そのもの(`generic_add_toFF`)。
2. **`n • ` と可換**なので `τ_T(n·生成点) = n·(生成点 + T) = n·生成点 + n·T = n·生成点`
   (`autFF_nsmul_generic`、`n·T = 0`)。
3. したがって **`τ(x_n) = x_n`、`τ(y_n) = y_n`**(`Point.some` の単射性)。
4. **環準同型は生成元で決まる**(第 119 の `coordinateRing_hom_ext`)ので
   `τ ∘ μ = μ`(`aut_comp_mulByN`)。

★★**`Point.map` が `→+` として mathlib にあった**ことが効いた——
群法則との可換性を自前で示す必要はない。

## ★★★残るのは素点の輸送だけ

`τ ∘ μ = μ` から `w_v(τ(μ f_P)) = w_v(μ f_P)`。★あとは

    w_v ∘ τ = w_{v'}      (Q_{v'} = Q_v + T)

を示せば `count_{v'}(μ f_P) = count_v(μ f_P)`、すなわち `e_v` の一定性が出る。
★★中心の同定は第 167 の `pullbackPrime_curveHom` がそのまま使える。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `generic_add_toFF` | ★★★★★平行移動した座標 = 生成点 + 定数点 |
| `autFF` / `autFF_some` | ★★★★★★自己同型が誘導する点の写像(`→+`) |
| `autFF_generic` | ★★★★★★**`τ(生成点) = 生成点 + T`** |
| `autFF_nsmul_generic` | ★★★★★★★**`n·T = 0` なら `n·生成点` を固定** |
| `aut_comp_mulByN` | ★★★★★★★**`τ ∘ μ = μ`** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine Polynomial IsDedekindDomain nonZeroDivisors

variable {F : Type} [Field F] [DecidableEq F] (W : WeierstrassCurve.Affine F) [W.IsElliptic]

/-! ## ★★★★★平行移動した座標は「生成点 + 定数点」 -/

/-- ★★★★★**平行移動した座標は「生成点 + 定数点」である**。

★`Point.add_some` そのもの。★★`slope` の `DecidableEq` 実体がずれるので、
傾きを `slopeFF` に書き換えてから合わせる。 -/
theorem generic_add_toFF {x₀ y₀ : F} (hQ : W.Nonsingular x₀ y₀) :
    ABC3.Found.GaloisRep.genericPoint W + toFF W (Point.some x₀ y₀ hQ)
      = Point.some (translateX W x₀ y₀) (translateY W x₀ y₀)
        (nonsingular_translate W hQ) := by
  rw [ABC3.Found.GaloisRep.genericPoint, toFF_some]
  have hne : ¬ (coordX W = algebraMap F W.FunctionField x₀
      ∧ coordY W = (W.map (algebraMap F W.FunctionField)).negY
          (algebraMap F W.FunctionField x₀) (algebraMap F W.FunctionField y₀)) :=
    fun hh => coordX_ne_const W x₀ hh.1
  have hs : (W.map (algebraMap F W.FunctionField)).slope (coordX W)
      (algebraMap F W.FunctionField x₀) (coordY W) (algebraMap F W.FunctionField y₀)
      = slopeFF W x₀ y₀ := by
    rw [WeierstrassCurve.Affine.slope_of_X_ne (coordX_ne_const W x₀)]; rfl
  rw [Point.add_some hne]
  simp only [hs, translateX, translateY]

/-! ## ★★★★★★自己同型が誘導する点の写像 -/

/-- ★★★★★★**関数体の `F`-自己同型が誘導する点の写像**(加法群の準同型)。 -/
noncomputable def autFF (τ : W.FunctionField ≃ₐ[F] W.FunctionField) :
    (W.map (algebraMap F W.FunctionField)).Point →+
      (W.map (algebraMap F W.FunctionField)).Point :=
  WeierstrassCurve.Affine.Point.map (W' := W) τ.toAlgHom

theorem autFF_some (τ : W.FunctionField ≃ₐ[F] W.FunctionField) {x y : W.FunctionField}
    (h : (W.map (algebraMap F W.FunctionField)).Nonsingular x y) :
    ∃ h' : (W.map (algebraMap F W.FunctionField)).Nonsingular (τ x) (τ y),
      autFF W τ (Point.some x y h) = Point.some (τ x) (τ y) h' :=
  ⟨_, WeierstrassCurve.Affine.Point.map_some _ _⟩

/-- ★★★★★★**平行移動の自己同型は生成点を `生成点 + T` に送る**。 -/
theorem autFF_generic (τ : W.FunctionField ≃ₐ[F] W.FunctionField)
    {x₀ y₀ : F} (hQ : W.Nonsingular x₀ y₀)
    (hx : τ (coordX W) = translateX W x₀ y₀) (hy : τ (coordY W) = translateY W x₀ y₀) :
    autFF W τ (ABC3.Found.GaloisRep.genericPoint W)
      = ABC3.Found.GaloisRep.genericPoint W + toFF W (Point.some x₀ y₀ hQ) := by
  obtain ⟨h', heq⟩ := autFF_some W τ (nonsingular_coord W)
  rw [generic_add_toFF W hQ, ABC3.Found.GaloisRep.genericPoint, heq]
  exact point_some_congr hx hy h' (nonsingular_translate W hQ)

/-- ★★★★★★★**`n·T = 0` なら自己同型は `n·生成点` を固定する**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★`Point.map` が `→+` なので `n • ` と可換。★★これが `[n] ∘ τ_T = [n]` の中身である。 -/
theorem autFF_nsmul_generic (τ : W.FunctionField ≃ₐ[F] W.FunctionField)
    {x₀ y₀ : F} (hQ : W.Nonsingular x₀ y₀)
    (hx : τ (coordX W) = translateX W x₀ y₀) (hy : τ (coordY W) = translateY W x₀ y₀)
    (n : ℕ) (hT : n • Point.some x₀ y₀ hQ = 0) :
    autFF W τ (n • ABC3.Found.GaloisRep.genericPoint W)
      = n • ABC3.Found.GaloisRep.genericPoint W := by
  rw [map_nsmul, autFF_generic W τ hQ hx hy, nsmul_add, ← map_nsmul (toFF W), hT, map_zero,
    add_zero]

/-- ★★★★★★★**`n·T = 0` なら `τ ∘ μ = μ`**——`μ f_P` は平行移動で不変。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★`τ` が `n·生成点` を固定することと、第 119 の `coordinateRing_hom_ext` から。 -/
theorem aut_comp_mulByN (τ : W.FunctionField ≃ₐ[F] W.FunctionField)
    {x₀ y₀ : F} (hQ : W.Nonsingular x₀ y₀)
    (hx : τ (coordX W) = translateX W x₀ y₀) (hy : τ (coordY W) = translateY W x₀ y₀)
    (n : ℕ) (hT : n • Point.some x₀ y₀ hQ = 0)
    (μ : W.CoordinateRing →+* W.FunctionField)
    (hμF : ∀ d : F, μ (algebraMap F W.CoordinateRing d) = algebraMap F W.FunctionField d)
    {xn yn : W.FunctionField}
    (hns : (W.map (algebraMap F W.FunctionField)).Nonsingular xn yn)
    (hμx : μ (genX W) = xn) (hμy : μ (genY W) = yn)
    (hμP : n • ABC3.Found.GaloisRep.genericPoint W = Point.some xn yn hns) :
    (τ.toAlgHom.toRingHom).comp μ = μ := by
  have hfix := autFF_nsmul_generic W τ hQ hx hy n hT
  rw [hμP] at hfix
  obtain ⟨h', heq⟩ := autFF_some W τ hns
  rw [heq] at hfix
  have hxe : τ xn = xn := by injection hfix
  have hye : τ yn = yn := by injection hfix
  refine coordinateRing_hom_ext _ _ (fun d => ?_) ?_ ?_
  · simp only [RingHom.comp_apply, hμF d, AlgHom.toRingHom_eq_coe, RingHom.coe_coe]
    exact τ.commutes d
  · simp only [RingHom.comp_apply, hμx, AlgHom.toRingHom_eq_coe, RingHom.coe_coe]
    exact hxe
  · simp only [RingHom.comp_apply, hμy, AlgHom.toRingHom_eq_coe, RingHom.coe_coe]
    exact hye

/-! ## ★出典の紐付け(`.src`) -/

def autFF_nsmul_generic.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の構成——平行移動が n·生成点を固定すること)",
    sectionId := "genell-thm-3-8" }

def aut_comp_mulByN.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の構成——μ f_P が平行移動で不変であること)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
