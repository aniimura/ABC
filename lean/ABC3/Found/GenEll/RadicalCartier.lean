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

★★★**茎で定義した有効 Cartier 因子(`IsEffectiveCartierStalk`)なら大域化が届いた。**
必要なのは「茎を取る操作と根基を取る操作が交換する」ことだけであり、
それは**茎が局所化である**ことから出る:

- `IsAffineOpen.isLocalization_stalk` —— アフィン開集合の切断から茎は局所化
- `IsLocalization.map_radical` —— 局所化は根基と交換する

★どちらも mathlib にある(2026-08-17 実測)。

## ★★これで (ii) はどこまで来たか

| 段 | 状態 |
|---|---|
| UFD の局所的な主張 | ★`RadicalPrincipal.lean` で取得済 |
| **茎と根基の交換** | ★★★**本ファイルで取得**(`stalk_radical`) |
| **茎の水準への持ち上げ** | ★★**本ファイルで取得**(`IsEffectiveCartierAt.radical`) |
| **scheme への大域化** | ★★★**本ファイルで取得**(`isEffectiveCartierStalk_radical`) |
| **被約性の定義と冪等性** | ★★**本ファイルで取得** |
| 正則局所環は UFD(Auslander–Buchsbaum) | ★★★**mathlib に無い** |

★★★**残るのは Auslander–Buchsbaum 1 本だけになった。**

本ファイルは正則性の代わりに **`UniqueFactorizationMonoid` を仮定**する——
正則 ⟹ UFD なので、これは原文より**弱い仮定で強い主張**である。
★ただし原文の「正則部分」という形で使うには A–B が要る。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory

/-! ## ★★茎と根基は交換する -/

/-- ★★★**摩擦を抜けた形** —— 点を **`↥U` の変数として持つ**。

## ★★★器具の記録 —— 点の持ち方だけで通ったり通らなかったりする

最初は `x : X` と `hxU : x ∈ U` から `⟨x, hxU⟩` を組んで
`hU.isLocalization_stalk ⟨x, hxU⟩` を使おうとした。★これは**通らない**——
その instance は `Algebra Γ(X,U) (stalk ↑⟨x, hxU⟩)` に付いており、
目標は `stalk x` に付いたものを要求する。
`↑⟨x, hxU⟩` と `x` は定義的に等しいのに、`IsLocalization` の型クラス引数を通すと
★**`whnf` が 100 万ヒートビートでも終わらない**(`erw` も `convert` も効かない)。

★★★**`y : ↥U` を変数として持つと通る。** instance が目標と**構文的に**合うからである。
★★数学は 1 ミリも変わっていない——**点の持ち方だけ**が違う。 -/
theorem ideal_radical_map_germ {X : Scheme} (I : X.IdealSheafData)
    {U : X.Opens} (hU : IsAffineOpen U) (y : U) :
    ((I.ideal ⟨U, hU⟩).radical).map (X.presheaf.germ U y.1 y.2).hom
      = ((I.ideal ⟨U, hU⟩).map (X.presheaf.germ U y.1 y.2).hom).radical := by
  haveI := hU.isLocalization_stalk y
  exact IsLocalization.map_radical (S := X.presheaf.stalk y.1)
    (hU.primeIdealOf y).asIdeal.primeCompl _

/-- ★★★**イデアル層の茎は根基と交換する**。

原文 (GenEll p.8):
> is also an effective Cartier divisor. We shall say that E is reduced if E = Ered.

★機構は「茎は局所化である」こと
(`IsAffineOpen.isLocalization_stalk` + `IsLocalization.map_radical`)。
★★点を `↥U` の変数として持つ形(`ideal_radical_map_germ`)を経由すれば通る。 -/
theorem stalk_radical {X : Scheme} (I : X.IdealSheafData) (x : X) :
    (Scheme.IdealSheafData.radical I).stalk x = (I.stalk x).radical := by
  obtain ⟨_, ⟨U, hU, rfl⟩, hxU, -⟩ :=
    X.isBasis_affineOpens.exists_subset_of_mem_open (Set.mem_univ x) isOpen_univ
  rw [Scheme.IdealSheafData.stalk_eq_map (U := ⟨U, hU⟩) hxU,
    Scheme.IdealSheafData.stalk_eq_map (U := ⟨U, hU⟩) hxU]
  exact ideal_radical_map_germ I hU ⟨x, hxU⟩

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

    (h : IsEffectiveCartierStalk I) (hne : ∀ x : X, I.stalk x ≠ ⊥) :
    IsEffectiveCartierStalk (Scheme.IdealSheafData.radical I) := by
  intro x
  rw [stalk_radical]
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

def stalk_radical.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 8,
    item := "Definition 1.5, (ii)(茎と根基の交換のみ)",
    sectionId := "genell-def-1-5" }

end ABC3.Found.GenEll
