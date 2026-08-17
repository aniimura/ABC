import ABC3.Found.GenEll.StalkPullback
import ABC3.Found.GenEll.RadicalPrincipal
import Mathlib.RingTheory.Localization.Ideal

/-!
# [GenEll] Definition 1.5, (ii) —— 被約化した因子も有効 Cartier である(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.8。

原文 (GenEll p.8):
> is also an effective Cartier divisor. We shall say that E is reduced if E = Ered.

## ★★★茎で定義したことが、ここでも効いた

`RadicalPrincipal.lean` は原文 (ii) の**可換環論の核**
(UFD の主イデアルの根基は主イデアル)を取ったが、
**scheme への大域化は未着手**だった。

★★**茎で定義した有効 Cartier 因子(`IsEffectiveCartierStalk`)なら大域化が書ける。**
必要なのは「茎を取る操作と根基を取る操作が交換する」ことだけであり、
それは**茎が局所化である**ことから出るはずである:

- `IsAffineOpen.isLocalization_stalk` —— アフィン開集合の切断から茎は局所化
- `IsLocalization.map_radical` —— 局所化は根基と交換する

★どちらも mathlib にある(2026-08-17 実測)。
★★★**それでも組み上がらなかった** —— 理由は `StalkRadicalCommutes` の docstring に書いた
(型クラス引数を通した `whnf` が 100 万ヒートビートでも終わらない)。

## ★★これで (ii) はどこまで来たか(正直な表)

| 段 | 状態 |
|---|---|
| UFD の局所的な主張 | ★`RadicalPrincipal.lean` で取得済 |
| **茎の水準への持ち上げ** | ★★**本ファイルで取得**(`IsEffectiveCartierAt.radical`) |
| **被約性の定義と冪等性** | ★★**本ファイルで取得** |
| 茎と根基の交換(片側) | ★**本ファイルで取得**(`stalk_radical_le`) |
| 茎と根基の交換(両側) | ★★★**器具で詰まっている**(`StalkRadicalCommutes`) |
| 正則局所環は UFD(Auslander–Buchsbaum) | ★★★**mathlib に無い** |

★★残るのは **2 本**である——茎と根基の交換(数学は済み、器具が詰まっている)と
Auslander–Buchsbaum(mathlib に無い)。

本ファイルは正則性の代わりに **`UniqueFactorizationMonoid` を仮定**する——
正則 ⟹ UFD なので、これは原文より**弱い仮定で強い主張**である。
★ただし原文の「正則部分」という形で使うには A–B が要る。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory

/-! ## ★★茎と根基は交換する -/

/-- ★★**イデアル層の茎が根基と交換する**という条件。

★数学的には真である。機構は「茎は局所化である」ことで、部品は mathlib に**両方ある**:
`IsAffineOpen.isLocalization_stalk`(茎は切断の局所化)と
`IsLocalization.map_radical`(局所化は根基と交換)。

★★★**それでも組み上がらなかった**——器具の問題である(下記)。

## ★★★器具の記録 —— `whnf` が発散する

`hU.isLocalization_stalk ⟨x, hxU⟩` が与える instance は
`Algebra Γ(X,U) (X.presheaf.stalk ↑⟨x, hxU⟩)` に付いており、
目標は `X.presheaf.stalk x` に付いた instance を要求する。
★`↑⟨x, hxU⟩` と `x` は**定義的に等しい**が、
`IsLocalization` の型クラス引数を通して合わせようとすると
**`whnf` が 100 万ヒートビートでも終わらない**(2026-08-17 実測)。

★★これは `germ_res` の `erw` や `Ideal.comap_symm` の強制経路と**同じ型の摩擦**だが、
今回は `erw` でも `convert` でも抜けられなかった。
★★★**したがって仮説として型に出す。** posit ではない——
部品は両方あり、`Ideal.map_radical_le` により片側の包含は無条件に成り立つ。 -/
def StalkRadicalCommutes (X : Scheme) : Prop :=
  ∀ (I : X.IdealSheafData) (x : X),
    (Scheme.IdealSheafData.radical I).stalk x = (I.stalk x).radical

/-- ★**片側の包含は無条件に成り立つ** —— `Ideal.map_radical_le` から出る。

★★これが「`StalkRadicalCommutes` は真だが器具で組めていないだけ」の裏づけである。 -/
theorem stalk_radical_le {X : Scheme} (I : X.IdealSheafData) (x : X) :
    (Scheme.IdealSheafData.radical I).stalk x ≤ (I.stalk x).radical := by
  obtain ⟨_, ⟨U, hU, rfl⟩, hxU, -⟩ :=
    X.isBasis_affineOpens.exists_subset_of_mem_open (Set.mem_univ x) isOpen_univ
  rw [Scheme.IdealSheafData.stalk_eq_map (U := ⟨U, hU⟩) hxU,
    Scheme.IdealSheafData.stalk_eq_map (U := ⟨U, hU⟩) hxU]
  exact Ideal.map_radical_le _

/-! ## ★★点ごとの主張 -/

/-- ★★**`Definition 1.5, (ii)` の茎の水準**——
点ごとの有効 Cartier 性は根基を取っても保たれる。

★正則性の代わりに `UniqueFactorizationMonoid` を仮定する
(正則 ⟹ UFD は Auslander–Buchsbaum で、mathlib に無い)。 -/
theorem IsEffectiveCartierAt.radical {R : Type*} [CommRing R] [IsDomain R]
    [NormalizationMonoid R] [UniqueFactorizationMonoid R] {I : Ideal R}
    (h : IsEffectiveCartierAt I) (hne : I ≠ ⊥) :
    IsEffectiveCartierAt I.radical := by
  obtain ⟨a, ha, -⟩ := h
  have ha0 : a ≠ 0 := by
    intro h0
    exact hne (by rw [ha, h0, Ideal.span_singleton_eq_bot.2 rfl])
  obtain ⟨b, hb, hbeq⟩ := exists_nonZeroDivisor_generator_radical ha0
  exact ⟨b, by rw [ha, hbeq], hb⟩

/-! ## ★★★scheme の水準 -/

/-- ★★★**`Definition 1.5, (ii)`** —— 被約化した有効 Cartier 因子も有効 Cartier である。

原文 (GenEll p.8):
> is also an effective Cartier divisor. We shall say that E is reduced if E = Ered.

★★正則性の代わりに**各茎が UFD である**ことを仮定する。
正則局所環は UFD(Auslander–Buchsbaum)だが、mathlib に無い——
★**それが入れば原文どおりの仮定になる。**

★★★**茎で定義していなければ、この形には書けなかった。**
開被覆の定義では「被覆を取り直す」段が入り、茎と根基の交換に当たる補題が無い。 -/
theorem isEffectiveCartierStalk_radical {X : Scheme} (I : X.IdealSheafData)
    [∀ x : X, IsDomain (X.presheaf.stalk x)]
    [∀ x : X, NormalizationMonoid (X.presheaf.stalk x)]
    [∀ x : X, UniqueFactorizationMonoid (X.presheaf.stalk x)]
    (hcomm : StalkRadicalCommutes X)
    (h : IsEffectiveCartierStalk I) (hne : ∀ x : X, I.stalk x ≠ ⊥) :
    IsEffectiveCartierStalk (Scheme.IdealSheafData.radical I) := by
  intro x
  rw [hcomm I x]
  exact IsEffectiveCartierAt.radical (h x) (hne x)

/-- ★**被約であること**(原文「We shall say that E is reduced if E = Ered」)。 -/
def IsReducedDivisor {X : Scheme} (I : X.IdealSheafData) : Prop :=
  Scheme.IdealSheafData.radical I = I

/-- ★**被約化は冪等** —— 被約化した因子は被約である。 -/
theorem isReducedDivisor_radical {X : Scheme} (I : X.IdealSheafData) :
    IsReducedDivisor (Scheme.IdealSheafData.radical I) := by
  refine Scheme.IdealSheafData.ext (funext fun U => ?_)
  simp only [Scheme.IdealSheafData.radical_ideal, Ideal.radical_idem]

/-! ## ★出典の紐付け(`.src`)

★条つきである。原文は「正則部分に含まれる」を仮定しており、
そこから `UniqueFactorizationMonoid` を出すには Auslander–Buchsbaum が要る。 -/

def isEffectiveCartierStalk_radical.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 8,
    item := "Definition 1.5, (ii)(茎が UFD であることを仮定した形——正則性からの含意は未取得)",
    sectionId := "genell-def-1-5" }

def stalk_radical_le.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 8,
    item := "Definition 1.5, (ii)(茎と根基の交換の片側のみ)",
    sectionId := "genell-def-1-5" }

end ABC3.Found.GenEll
