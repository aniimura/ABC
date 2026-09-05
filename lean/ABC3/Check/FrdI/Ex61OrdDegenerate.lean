/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Skeleton.Divisor.SchemeWeil

/-!
# [FrdI] Example 6.1 の因子論——現行 statement は**零写像で充足される**

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.109(`Example 6.1`)。

原文 (FrdI p.109):
> normal [geometrically integral] variety over a field k; K the function field of V ;

原文 (FrdI p.109):
> DK a set of Q-Cartier prime divisors on V . The connected objects of the Galois category

## ★★なぜこの検査が要るか

`Skeleton/Divisor/SchemeWeil.lean` の `ordAtDiv` は **`sorry` 本体の `def`** である。
`def` が `sorry` だということは、その中身について**何も分かっていない**ということで、
`ordAtDiv` を縛るのは同じファイルの 3 定理

| 定理 | 内容 |
|---|---|
| `ordAtDiv_mul` | `ord_x(fg) = ord_x(f) + ord_x(g)`(`f ≠ 0`, `g ≠ 0`) |
| `finite_support_ordAtDiv` | `{x | ord_x(f) ≠ 0}` は有限(`f ≠ 0`) |
| `divOfFn_mul` | `div(fg) = div(f) + div(g)` |

だけである。★**この 3 条は `ord ≡ 0`・`div ≡ 0` で全部成り立つ**
(`0 = 0 + 0`、`∅.Finite`、`0 = 0 + 0`)。

したがって `SchemeWeil.lean` の 6 つの `sorry` のうち

* `ordAtDiv`(`def`)・`ordAtDiv_mul`・`finite_support_ordAtDiv`・
  `divOfFn`(`def`)・`divOfFn_mul` の **5 つは数学的内容ゼロで「正直に」埋まる**
* 本物の数学が要るのは `isDiscreteValuationRing_stalk_of_codimOne` の **1 件だけ**

★★これは「`Skeleton` に嘘がある」という話ではない。**現行の statement 集合が
`ord` を一意に決めていない**という話である。本ファイルはその退化を Lean で
**構成して固定する**(`zeroWeilOrd`)。

## ★★`hnorm : IsNormalScheme X` を足すだけでは退化は消えない

`ordAtDiv` 以降の 5 つからは `IsNormalScheme X` が**丸ごと抜けている**
(`isDiscreteValuationRing_stalk_of_codimOne` にしか出てこない)。
これは一見「正規性を足せば直る」ように見えるが、**直らない**:

`WeilOrdSpecOldNormal`(3 条すべてに `hnorm : IsNormalScheme X` を足した形)も
やはり零写像で充足される(`example_6_1_ord_satisfied_by_zero_even_with_normal`)。

★理屈はこう。`hnorm` は「余次元 1 の茎が DVR になり、**正直な `ord` が定義できる**」
ための条件であって、「**零写像を排除する**」ための条件ではない。
前者は存在の話、後者は一意性の話で、片方をいくら強めても他方は出ない。

## ★★排除には「錨」が要る

退化を殺すのは仮定ではなく**値を 1 点で固定する条件**である:

```
HasOrdAnchor S : ∀ X [IsIntegral X] (x : PrimeDivisorPt X), ∃ f, S.ord X x f = 1
```

これは零写像を即座に排除する——正確には、零写像がこれを満たすのは
**素因子が 1 つも無いとき**に限る(`zeroWeilOrd_not_hasOrdAnchor`)。

★並行して別の agent が `Found/Divisor/SchemeWeilOrd.lean` に
`exists_ordPt_eq_one`(素元 `π` を取れば `ordPt = 1`)を立てている。
それが入れば、`Skeleton` 側は `ordAtDiv` に錨を 1 本足すだけで退化が閉じる。

## ★★下流への伝播(`Skeleton/Divisor/Cartier/Example61.lean`)

`IsCartierDiv` は `ordAtDiv` を**直に**使う:

```
IsCartierDiv X D := ∀ x : X, ∃ U ∋ x, ∃ f : K(X)ˣ, ∀ y ∈ U, D y = ordAtDiv X y f
```

`ord ≡ 0` だと右辺が消えるので、`IsCartierDiv X D ↔ D = 0` になり
(`isCartierDivOld_zeroWeilOrd_iff`)、`cartierSubgroup X = ⊥`
(`cartierSubgroupZero_eq_bot`)。したがって `Cartier/Example61.lean` の 4 つの
`sorry` も**すべて零写像で閉じる**:

| `Cartier/Example61.lean` の `sorry` | 退化での埋まり方 |
|---|---|
| `isCartierDiv_zero` | `isCartierDivOld_zero` |
| `isCartierDiv_add` | `isCartierDivOld_add` |
| `isCartierDiv_neg` | `isCartierDivOld_neg` |
| `isQCartierSubgroup_of_forall_isQCartier` | ★前提 `hQ` が**矛盾**する(`isQCartier_hypothesis_absurd_under_zeroWeilOrd`)ので空虚に真 |

★最後の行が効いている。`IsQCartierDiv X (single (ι s) 1)` は退化の下では
`single (ι s) 1 = 0`、すなわち `(1 : ℤ) = 0` を意味する。つまり
**`Q`-Cartier 性の仮定そのものが退化を排除する錨になっている**——ただし
`Cartier/Example61.lean` はその仮定を `isQCartierSubgroup_of_forall_isQCartier`
の中でしか使わないので、他の 3 件は救われない。

## ★★さらに下流(`Skeleton/Divisor/Cartier/Theorem62.lean`)

`pullbackCartier` も `sorry` 本体の `def` で、縛るのは `pullbackCartier_add`・
`isCartierDiv_pullbackCartier`・`pullbackCartier_nonneg` の 3 つだけ。
★こちらは `ord` の退化すら要らず、**`pullbackCartier ≡ 0` と置くだけ**で 4 つとも通る
(`0 = 0 + 0`・`IsCartierDiv Y 0`・`0 ≤ 0`)。`theorem_6_2_pullback_satisfied_by_zero`。

★★合計 **5 + 4 + 4 = 13 件**が零写像だけで閉じることを、本ファイルで実際に構成した
(`Skeleton/Divisor/` の `weil` / `cartier` 両鎖の `sorry` は 14 件で、残る 1 件が
本物の数学を要する `isDiscreteValuationRing_stalk_of_codimOne` である)。

## ★この検査の位置づけ

「落とした条件は、主張を偽にするか自明にするかのどちらかになる」例の **8 つ目**
(`InertiaDegeneracy`・`Theorem42Degenerate`・`Def32Degenerate`・`Cor33Degenerate`・
`Prop22Degenerate`・`Cor13Degenerate`・本件。`VLocFalse` は枠そのものの反証なので
この通番には入れていない)。

★本件は `InertiaDegeneracy`(`I_K := ⊥` で `Corollary 1.3` が通る)と**同型の病**である
——`Prop22Degenerate` / `Cor13Degenerate` のように「**偽**」を示すのではなく、
「**自明**」を示す。`lake build` も `sorry` 検査も `#print axioms` も通ってしまう型で、
機械では捕まらない。だから**こうして退化 witness を明示的に構成して固定する**。
-/

namespace ABC3.Check.FrdI

open AlgebraicGeometry CategoryTheory ABC3.Skeleton.Divisor

universe u

/-! ## ★1. 旧形——`Skeleton/Divisor/SchemeWeil.lean` が `ord`・`div` に課している条件の全部

★`Skeleton` の `ordAtDiv` / `divOfFn` は `sorry` 本体の `def` なので、
それ自身については何も証明できない。そこで**課されている条件だけを構造体に括り出し**、
「その条件を満たす組」を走らせる形に書き直す(`Prop22Degenerate` の
`RecoverableAsAddModuleOld` と同じ手口)。 -/

/-- **旧形** —— `SchemeWeil.lean` の `ordAtDiv` / `divOfFn` と、それらを縛る 3 定理。
★★フィールドはこれで**全部**である。 -/
structure WeilOrdSpecOld where
  /-- `ord_x : K(X) → ℤ`(`Skeleton.Divisor.ordAtDiv`)。 -/
  ord : ∀ (X : Scheme.{u}) [IsIntegral X], PrimeDivisorPt X → X.functionField → ℤ
  /-- `div : K(X)^× → WeilDiv X`(`Skeleton.Divisor.divOfFn`)。 -/
  div : ∀ (X : Scheme.{u}) [IsIntegral X] [AlgebraicGeometry.IsNoetherian X],
    (X.functionField)ˣ → WeilDiv X
  /-- `Skeleton.Divisor.ordAtDiv_mul`。 -/
  ord_mul : ∀ (X : Scheme.{u}) [IsIntegral X] (x : PrimeDivisorPt X)
    (f g : X.functionField), f ≠ 0 → g ≠ 0 →
      ord X x (f * g) = ord X x f + ord X x g
  /-- `Skeleton.Divisor.finite_support_ordAtDiv`。 -/
  finite_support : ∀ (X : Scheme.{u}) [IsIntegral X] [AlgebraicGeometry.IsNoetherian X]
    (f : X.functionField), f ≠ 0 → {x : PrimeDivisorPt X | ord X x f ≠ 0}.Finite
  /-- `Skeleton.Divisor.divOfFn_mul`。 -/
  div_mul : ∀ (X : Scheme.{u}) [IsIntegral X] [AlgebraicGeometry.IsNoetherian X]
    (f g : (X.functionField)ˣ), div X (f * g) = div X f + div X g

/-- ★★**零写像は旧形を充足する** —— `0 = 0 + 0`・`∅.Finite`・`0 = 0 + 0`。 -/
def zeroWeilOrd : WeilOrdSpecOld.{u} where
  ord _ _ _ _ := 0
  div _ _ _ _ := 0
  ord_mul := by intros; simp
  finite_support := by intro X _ _ f _; simp
  div_mul := by intros; simp

@[simp] theorem zeroWeilOrd_ord (X : Scheme.{u}) [IsIntegral X] (x : PrimeDivisorPt X)
    (f : X.functionField) : zeroWeilOrd.ord X x f = 0 := rfl

@[simp] theorem zeroWeilOrd_div (X : Scheme.{u}) [IsIntegral X]
    [AlgebraicGeometry.IsNoetherian X] (f : (X.functionField)ˣ) :
    zeroWeilOrd.div X f = 0 := rfl

/-- ★★★★★★**[FrdI] Example 6.1 の現行 statement は零写像で充足される**。

`SchemeWeil.lean` の 5 つの `sorry`(`ordAtDiv`・`ordAtDiv_mul`・
`finite_support_ordAtDiv`・`divOfFn`・`divOfFn_mul`)は、
**数学的内容ゼロでいっぺんに埋まる**。 -/
theorem example_6_1_ord_satisfied_by_zero :
    ∃ S : WeilOrdSpecOld.{u},
      (∀ (X : Scheme.{u}) [IsIntegral X] (x : PrimeDivisorPt X) (f : X.functionField),
        S.ord X x f = 0) ∧
      (∀ (X : Scheme.{u}) [IsIntegral X] [AlgebraicGeometry.IsNoetherian X]
        (f : (X.functionField)ˣ), S.div X f = 0) := by
  refine ⟨zeroWeilOrd, ?_, ?_⟩ <;> intros <;> rfl

/-! ## ★2. `IsNormalScheme` を足しても退化は消えない

★`ordAtDiv` 以降の 5 つからは `IsNormalScheme X` が抜けている。
では足せば直るか——**直らない**。 -/

/-- **旧形 + 正規性** —— 3 条すべてに `hnorm : IsNormalScheme X` を足した形。 -/
structure WeilOrdSpecOldNormal where
  /-- `ord_x : K(X) → ℤ`(正規性つき)。 -/
  ord : ∀ (X : Scheme.{u}) [IsIntegral X], IsNormalScheme X → PrimeDivisorPt X →
    X.functionField → ℤ
  /-- `div : K(X)^× → WeilDiv X`(正規性つき)。 -/
  div : ∀ (X : Scheme.{u}) [IsIntegral X] [AlgebraicGeometry.IsNoetherian X],
    IsNormalScheme X → (X.functionField)ˣ → WeilDiv X
  /-- `ord` の乗法性(正規性つき)。 -/
  ord_mul : ∀ (X : Scheme.{u}) [IsIntegral X] (hnorm : IsNormalScheme X)
    (x : PrimeDivisorPt X) (f g : X.functionField), f ≠ 0 → g ≠ 0 →
      ord X hnorm x (f * g) = ord X hnorm x f + ord X hnorm x g
  /-- 台の有限性(正規性つき)。 -/
  finite_support : ∀ (X : Scheme.{u}) [IsIntegral X] [AlgebraicGeometry.IsNoetherian X]
    (hnorm : IsNormalScheme X) (f : X.functionField), f ≠ 0 →
      {x : PrimeDivisorPt X | ord X hnorm x f ≠ 0}.Finite
  /-- `div` の準同型性(正規性つき)。 -/
  div_mul : ∀ (X : Scheme.{u}) [IsIntegral X] [AlgebraicGeometry.IsNoetherian X]
    (hnorm : IsNormalScheme X) (f g : (X.functionField)ˣ),
      div X hnorm (f * g) = div X hnorm f + div X hnorm g

/-- ★正規性を足しても零写像は生き残る。 -/
def zeroWeilOrdNormal : WeilOrdSpecOldNormal.{u} where
  ord _ _ _ _ _ := 0
  div _ _ _ _ _ := 0
  ord_mul := by intros; simp
  finite_support := by intro X _ _ _ f _; simp
  div_mul := by intros; simp

/-- ★★★★**`hnorm : IsNormalScheme X` を足しても退化は消えない**。

★`hnorm` は「正直な `ord` が**定義できる**」ための条件(存在)であって、
「零写像を**排除する**」ための条件(一意性)ではない。 -/
theorem example_6_1_ord_satisfied_by_zero_even_with_normal :
    ∃ S : WeilOrdSpecOldNormal.{u},
      ∀ (X : Scheme.{u}) [IsIntegral X] (hnorm : IsNormalScheme X) (x : PrimeDivisorPt X)
        (f : X.functionField), S.ord X hnorm x f = 0 :=
  ⟨zeroWeilOrdNormal, by intros; rfl⟩

/-! ## ★3. 錨——退化を殺すのはこれである -/

/-- ★**錨** —— 各素因子で値 `1` を取る有理函数がある。

★正直な `ord` はこれを満たす(素元 `π` を取ればよい)。
`Found/Divisor/SchemeWeilOrd.lean` の `ordPt` について、これを
`exists_ordPt_eq_one` として立てる作業が並行して走っている。 -/
def HasOrdAnchor (S : WeilOrdSpecOld.{u}) : Prop :=
  ∀ (X : Scheme.{u}) [IsIntegral X] (x : PrimeDivisorPt X), ∃ f : X.functionField, S.ord X x f = 1

/-- ★★**錨は零写像を排除する** —— 零写像が錨を満たすのは、素因子が
1 つも無いとき **に限る**。 -/
theorem zeroWeilOrd_not_hasOrdAnchor (h : HasOrdAnchor zeroWeilOrd.{u})
    (X : Scheme.{u}) [IsIntegral X] : IsEmpty (PrimeDivisorPt X) := by
  refine ⟨fun x => ?_⟩
  obtain ⟨f, hf⟩ := h X x
  rw [zeroWeilOrd_ord] at hf
  exact zero_ne_one hf

/-! ## ★4. 下流——`Skeleton/Divisor/Cartier/Example61.lean` も一緒に潰れる -/

/-- **旧形の Cartier 因子** —— `Skeleton.Divisor.IsCartierDiv` の `ordAtDiv` を
`S.ord` に差し替えたもの(定義は逐語で同じ)。 -/
def IsCartierDivOld (S : WeilOrdSpecOld.{u}) (X : Scheme.{u}) [IsIntegral X]
    (D : WeilDiv X) : Prop :=
  ∀ x : X, ∃ (U : X.Opens) (_ : x ∈ U) (f : (X.functionField)ˣ),
    ∀ y : PrimeDivisorPt X, y.1 ∈ U → D y = S.ord X y (f : X.functionField)

/-- ★★**退化すると Cartier 性は「`D = 0`」に潰れる**。

→ の向き: 素因子 `y` に対して `x := y.1` を取れば、その近傍で `D y = 0`。
← の向き: `U := ⊤`・`f := 1` を取ればよい。 -/
theorem isCartierDivOld_zeroWeilOrd_iff (X : Scheme.{u}) [IsIntegral X] (D : WeilDiv X) :
    IsCartierDivOld zeroWeilOrd X D ↔ D = 0 := by
  constructor
  · intro h
    ext y
    obtain ⟨U, hU, f, hf⟩ := h y.1
    simpa using hf y hU
  · rintro rfl x
    exact ⟨⊤, trivial, 1, fun y _ => by simp⟩

/-- `Cartier/Example61.lean::isCartierDiv_zero` の退化版。 -/
theorem isCartierDivOld_zero (X : Scheme.{u}) [IsIntegral X] :
    IsCartierDivOld zeroWeilOrd X 0 :=
  (isCartierDivOld_zeroWeilOrd_iff X 0).mpr rfl

/-- `Cartier/Example61.lean::isCartierDiv_add` の退化版。 -/
theorem isCartierDivOld_add (X : Scheme.{u}) [IsIntegral X] {D E : WeilDiv X}
    (hD : IsCartierDivOld zeroWeilOrd X D) (hE : IsCartierDivOld zeroWeilOrd X E) :
    IsCartierDivOld zeroWeilOrd X (D + E) := by
  rw [isCartierDivOld_zeroWeilOrd_iff] at hD hE ⊢
  rw [hD, hE, add_zero]

/-- `Cartier/Example61.lean::isCartierDiv_neg` の退化版。 -/
theorem isCartierDivOld_neg (X : Scheme.{u}) [IsIntegral X] {D : WeilDiv X}
    (hD : IsCartierDivOld zeroWeilOrd X D) : IsCartierDivOld zeroWeilOrd X (-D) := by
  rw [isCartierDivOld_zeroWeilOrd_iff] at hD ⊢
  rw [hD, neg_zero]

/-- `Cartier/Example61.lean::cartierSubgroup` の退化版。 -/
def cartierSubgroupZero (X : Scheme.{u}) [IsIntegral X] : AddSubgroup (WeilDiv X) where
  carrier := {D | IsCartierDivOld zeroWeilOrd X D}
  zero_mem' := isCartierDivOld_zero X
  add_mem' hD hE := isCartierDivOld_add X hD hE
  neg_mem' hD := isCartierDivOld_neg X hD

/-- ★★★★**`Φ(L)^gp` に相当する群が `⊥` に潰れる** —— `Example 6.1` の単系論は
退化の下では**中身が無い**。 -/
theorem cartierSubgroupZero_eq_bot (X : Scheme.{u}) [IsIntegral X] :
    cartierSubgroupZero X = ⊥ := by
  ext D
  simpa [cartierSubgroupZero] using isCartierDivOld_zeroWeilOrd_iff X D

/-- ★退化の下では `Q`-Cartier 性も「`D = 0`」に潰れる(`WeilDiv` は捩れが無い)。 -/
theorem isQCartierDivOld_zeroWeilOrd_iff (X : Scheme.{u}) [IsIntegral X] (D : WeilDiv X) :
    (∃ n : ℕ, 0 < n ∧ IsCartierDivOld zeroWeilOrd X (n • D)) ↔ D = 0 := by
  constructor
  · rintro ⟨n, hn, h⟩
    have hnD : n • D = 0 := (isCartierDivOld_zeroWeilOrd_iff X _).mp h
    ext y
    have h2 : (n : ℤ) * D y = 0 := by
      have := congrArg (fun F : WeilDiv X => F y) hnD
      simpa [nsmul_eq_mul] using this
    have hn' : (n : ℤ) ≠ 0 := Int.natCast_ne_zero.mpr hn.ne'
    simpa using (mul_eq_zero.mp h2).resolve_left hn'
  · rintro rfl
    exact ⟨1, one_pos, (isCartierDivOld_zeroWeilOrd_iff X _).mpr (by simp)⟩

/-- ★★★**`isQCartierSubgroup_of_forall_isQCartier` の前提は退化の下で矛盾する**。

`Cartier/Example61.lean` の `_hQ : ∀ s, IsQCartierDiv X (Finsupp.single (ι s) 1)` は、
`ord ≡ 0` の下では `Finsupp.single (ι s) 1 = 0`、すなわち `(1 : ℤ) = 0` を意味する。
★つまりその `sorry` も(空虚に)埋まる。★同時にこれは
「`Q`-Cartier 性の仮定それ自体が錨になっている」ことを示している——ただし
`Cartier/Example61.lean` はこの仮定を他の 3 定理では使わないので、
そちらは救われない。 -/
theorem isQCartier_hypothesis_absurd_under_zeroWeilOrd
    (X : Scheme.{u}) [IsIntegral X] {S : Type u} (ι : S → PrimeDivisorPt X) (s : S)
    (hQ : ∃ n : ℕ, 0 < n ∧
      IsCartierDivOld zeroWeilOrd X (n • Finsupp.single (ι s) (1 : ℤ))) : False := by
  have h := (isQCartierDivOld_zeroWeilOrd_iff X _).mp hQ
  rw [Finsupp.single_eq_zero] at h
  exact one_ne_zero h

/-! ## ★5. さらに下流——`Skeleton/Divisor/Cartier/Theorem62.lean` の引き戻し

★こちらは `ord` の退化すら要らない。**`pullbackCartier ≡ 0` と置くだけ**で
4 定理が全部通る(`0 = 0 + 0`・`IsCartierDiv Y 0`・`0 ≤ 0`)。 -/

/-- **旧形** —— `Theorem62.lean` の `pullbackCartier` と、それを縛る 3 定理。
★フィールドはこれで全部である。 -/
structure PullbackCartierSpecOld where
  /-- `Skeleton.Divisor.pullbackCartier`。 -/
  pull : ∀ (X Y : Scheme.{u}) [IsIntegral X] [AlgebraicGeometry.IsNoetherian X]
    [IsIntegral Y] [AlgebraicGeometry.IsNoetherian Y],
    (Y ⟶ X) → (D : WeilDiv X) → IsCartierDivOld zeroWeilOrd X D → WeilDiv Y
  /-- `Skeleton.Divisor.pullbackCartier_add`。 -/
  pull_add : ∀ (X Y : Scheme.{u}) [IsIntegral X] [AlgebraicGeometry.IsNoetherian X]
    [IsIntegral Y] [AlgebraicGeometry.IsNoetherian Y] (ψ : Y ⟶ X) (D E : WeilDiv X)
    (hD : IsCartierDivOld zeroWeilOrd X D) (hE : IsCartierDivOld zeroWeilOrd X E),
    pull X Y ψ (D + E) (isCartierDivOld_add X hD hE)
      = pull X Y ψ D hD + pull X Y ψ E hE
  /-- `Skeleton.Divisor.isCartierDiv_pullbackCartier`。 -/
  pull_cartier : ∀ (X Y : Scheme.{u}) [IsIntegral X] [AlgebraicGeometry.IsNoetherian X]
    [IsIntegral Y] [AlgebraicGeometry.IsNoetherian Y] (ψ : Y ⟶ X) (D : WeilDiv X)
    (hD : IsCartierDivOld zeroWeilOrd X D),
    IsCartierDivOld zeroWeilOrd Y (pull X Y ψ D hD)
  /-- `Skeleton.Divisor.pullbackCartier_nonneg`。 -/
  pull_nonneg : ∀ (X Y : Scheme.{u}) [IsIntegral X] [AlgebraicGeometry.IsNoetherian X]
    [IsIntegral Y] [AlgebraicGeometry.IsNoetherian Y] (ψ : Y ⟶ X) (D : WeilDiv X)
    (hD : IsCartierDivOld zeroWeilOrd X D), 0 ≤ D → 0 ≤ pull X Y ψ D hD

/-- ★★**零引き戻しは旧形を充足する**。 -/
def zeroPullbackCartier : PullbackCartierSpecOld.{u} where
  pull _ _ _ _ _ _ _ _ _ := 0
  pull_add := by intros; simp
  pull_cartier X Y _ _ _ _ _ _ _ := isCartierDivOld_zero Y
  pull_nonneg := by intros; exact le_refl 0

/-- ★★★★**`Theorem62.lean` の 4 つの `sorry` も零写像で閉じる**。

★ここは `ord` の退化とは独立で、`pullbackCartier ≡ 0` だけで通る。 -/
theorem theorem_6_2_pullback_satisfied_by_zero :
    ∃ P : PullbackCartierSpecOld.{u},
      ∀ (X Y : Scheme.{u}) [IsIntegral X] [AlgebraicGeometry.IsNoetherian X]
        [IsIntegral Y] [AlgebraicGeometry.IsNoetherian Y] (ψ : Y ⟶ X) (D : WeilDiv X)
        (hD : IsCartierDivOld zeroWeilOrd X D), P.pull X Y ψ D hD = 0 :=
  ⟨zeroPullbackCartier, by intros; rfl⟩

end ABC3.Check.FrdI
