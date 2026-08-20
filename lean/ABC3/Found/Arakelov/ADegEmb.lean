import ABC3.Found.Arakelov.APicWitness
import Mathlib.NumberTheory.NumberField.InfinitePlace.Embeddings

/-!
# Arakelov (D2) 第 300 ブロック —— **正規化次数は埋め込みの平均**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.4。

原文 (GenEll p.4):
> — where xF : Spec(OF ) → X is any morphism that gives rise to x.

## ★★★★★★`deg_F` は**純アルキメデス**でよい

`deg_F : APic(Spec 𝓞_F) → ℝ` は群準同型であり、
★`Pic(Spec 𝓞_F) ≅ Cl(F)` は**有限**なので、`ℝ` は捻れを持たず
**幾何部分の寄与は 0 に強制される**。
★★★したがって `deg_F` は Green 関数だけで書ける。

## ★★★素点でなく**埋め込み**で和を取る

素点 `v` で和を取ると、複素素点で 2 つの共役埋め込みのうち**一方だけ**を選ぶことになり、
★★底変換不変性が**共役不変性の仮定なしには落ちる**
(反例: `F = ℚ(i)`、`φ = 複素共役`)。

★★★**埋め込みで和を取れば**、共役の対がそのまま両方入るので

    deg_F(L̄) = -(1/[F:ℚ]) Σ_{σ : F →+* ℂ} g(p_σ)

は**無条件で**共役対称であり、底変換でも不変になる
(`Fintype.card (F →+* ℂ) = [F:ℚ]`、拡大では各 σ がちょうど `[K:F]` 回現れる)。

## ★★本ブロックで取れるもの

| 定義・定理 | 内容 |
|---|---|
| `embPoint` | ★埋め込みが定める複素点 |
| `degFOf` | ★★★正規化次数 |
| `degFOf_mul` | ★★加法性 |
| `degFOf_scale` | ★★★★**計量の定数倍で `-c` 動く**(退化を殺す欄) |
-/

namespace ABC3.Interface.Arakelov

open AlgebraicGeometry CategoryTheory NumberField ABC3.Found.Arakelov

attribute [local instance] aPicGroup

variable (F : Type) [Field F] [NumberField F]

/-- ★**埋め込み `σ : F →+* ℂ` が定める `Spec 𝓞_F` の複素点**。 -/
noncomputable def embPoint (σ : F →+* ℂ) :
    Spec (CommRingCat.of ℂ) ⟶ Spec (CommRingCat.of (𝓞 F)) :=
  Spec.map (CommRingCat.ofHom (σ.comp (algebraMap (𝓞 F) F)))

/-- ★★★**正規化次数**——Green 関数の、埋め込みにわたる平均の符号違い。 -/
noncomputable def degFOf
    (x : picardDataWitness.Pic (Spec (CommRingCat.of (𝓞 F)))
      × Multiplicative (arcCM (Spec (CommRingCat.of (𝓞 F))))) : ℝ :=
  -(∑ σ : F →+* ℂ, (x.2 : arcCM (Spec (CommRingCat.of (𝓞 F)))).fn (embPoint F σ))
    / (Module.finrank ℚ F : ℝ)

/-- ★★**加法性**。 -/
theorem degFOf_mul
    (x y : picardDataWitness.Pic (Spec (CommRingCat.of (𝓞 F)))
      × Multiplicative (arcCM (Spec (CommRingCat.of (𝓞 F))))) :
    degFOf F (x * y) = degFOf F x + degFOf F y := by
  show -(∑ σ : F →+* ℂ, ((x.2 * y.2 : Multiplicative _) :
      arcCM (Spec (CommRingCat.of (𝓞 F)))).fn (embPoint F σ)) / _ = _
  have hs : ∀ σ : F →+* ℂ,
      ((x.2 * y.2 : Multiplicative (arcCM (Spec (CommRingCat.of (𝓞 F))))) :
          arcCM (Spec (CommRingCat.of (𝓞 F)))).fn (embPoint F σ)
        = (x.2 : arcCM _).fn (embPoint F σ) + (y.2 : arcCM _).fn (embPoint F σ) :=
    fun _ => rfl
  simp only [hs, Finset.sum_add_distrib, degFOf]
  field_simp
  ring

/-- ★★★★**計量の定数倍で次数は `-c` 動く**。

原文 (GenEll p.4):
> — where xF : Spec(OF ) → X is any morphism that gives rise to x.

★★★これが「`deg_F ≡ 0`」の退化を殺す欄である。
★`Fintype.card (F →+* ℂ) = [F:ℚ]` がちょうど正規化の分母を消す。 -/
theorem degFOf_scale
    (L : picardDataWitness.Pic (Spec (CommRingCat.of (𝓞 F))))
    (m : TorsorMetric (Spec (CommRingCat.of (𝓞 F)))
      (picardDataWitness.sheafOf (Spec (CommRingCat.of (𝓞 F))) L)) (c : ℝ) :
    degFOf F (aPicDataWitness.ofMetric _ L (TorsorMetric.scale c m))
      = degFOf F (aPicDataWitness.ofMetric _ L m) - c := by
  have hcard : (Fintype.card (F →+* ℂ) : ℝ) = (Module.finrank ℚ F : ℝ) := by
    exact_mod_cast congrArg (Nat.cast (R := ℝ)) (Embeddings.card F ℂ)
  have hpos : (Module.finrank ℚ F : ℝ) ≠ 0 := by
    have := Module.finrank_pos (R := ℚ) (M := F)
    positivity
  show -(∑ σ : F →+* ℂ, (m.green (embPoint F σ) + c)) / (Module.finrank ℚ F : ℝ)
    = -(∑ σ : F →+* ℂ, m.green (embPoint F σ)) / (Module.finrank ℚ F : ℝ) - c
  rw [Finset.sum_add_distrib, Finset.sum_const, Finset.card_univ, nsmul_eq_mul, hcard]
  field_simp
  ring

/-! ## ★出典の紐付け(`.src`) -/

def degFOf_scale.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 4,
    item := "Definition 1.1, (ii)(正規化次数 deg_F が計量に依存すること)",
    sectionId := "genell-def-1-1-ii" }

end ABC3.Interface.Arakelov
