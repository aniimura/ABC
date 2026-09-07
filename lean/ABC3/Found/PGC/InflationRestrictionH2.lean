import ABC3.Found.PGC.GroupCohomologyFinite

/-!
# [pGC] Proposition 1.1 —— `H²` の inflation-restriction 完全列(mathlib 側の欠落)

`Found/PGC/GroupCohomologyFinite.lean` は Proposition 1.1 の離散側の欠落 (a)(b) を埋め、
残る (c)(連続コホモロジーと双対性)について **4 つ**を名指しした。本ファイルはその 1 つ目、

> 1. `H²` の inflation-restriction 完全列
>    (`H¹(G/S,A^S) → H¹(G,A) → H¹(S,A)^{G/S} → H²(G/S,A^S) → H²(G,A)`)

を埋める。★mathlib は `H1InfRes` / `H1InfRes_exact`(次数 1 だけ)しか持たない。

## ★★本ファイルの到達点(射程)

`S ⊴ G`、`A : Rep k G` について、短複体

  `H²(G ⧸ S, A^S) --inf--> H²(G, A) --res--> H²(S, A)`

を `H2InfRes` として定義し、**`H¹(S, A) = 0` のもとで**

1. `inf` が**単射**(`inflation₂_injective` / `mono_H2InfRes_f` / `mono_inflationH2`)
   ★これは 5 項完全列の**右端の完全性**(`H²(G/S,A^S)` における完全性)そのものである。
   `H¹(S,A) = 0` なら transgression の source が消えるので、右端の完全性は `inf` の単射性に一致する。
2. 短複体 `H2InfRes` が **`H²(G,A)` において完全**(`H2InfRes_exact`)

の 2 つを証明した。★**すなわち `0 → H²(G/S,A^S) → H²(G,A) → H²(S,A)` が完全**である
(古典的な「`H^i(S,A) = 0` for `i < n` ならば `0 → H^n(G/S,A^S) → H^n(G,A) → H^n(S,A)` が完全」
の `n = 2` の場合)。

★おまけに mathlib に無い `map₂_one`(`map (1 : G →* H) φ 2 = 0`)も埋めた
(mathlib は `map₁_one` しか持たない)。

## ★★仮定に置いたもの(名指し)

★**唯一の仮定は `H¹(S, A) = 0`** である。Lean では
`Subsingleton ↑(groupCohomology (Rep.res S.subtype A) 1)` と書く
(圏論的な `IsZero` 版も `..._of_isZero` として用意した)。

* ★**`axiom` も `structure` も `sorry` も置いていない。**
* ★**「通すための偽の仮説」は作っていない。** 仮定 `H¹(S,A) = 0` は
  5 項完全列の transgression の source を消すための**古典的に必要な仮定**であり、
  ★退化の自己検査 `inflation₂_bijective_bot`(`S = ⊥` のとき仮定は自動的に成り立ち、
  しかも `inf` は**全単射**になる)で空虚でないことを確認した。
* ★消費側(`Γ_K` の有限商の塔)がこの仮定を供給できるかは**別ノード**である。

## ★★★測定 —— mathlib に何が無いか(2026-09-07、コマンドと出力)

```
awk -F'\t' '$2 ~ /[Ii]nf(Res|NatTrans|lation)/ {print $2"\t"$3}' .cache/mathlib-index.txt
```
→ `groupCohomology.H1InfRes`(Functoriality.lean:342) /
  `groupCohomology.H1InfRes_exact`(:370) / `groupCohomology.infNatTrans`(:534) /
  `groupHomology.coinfNatTrans` の 4 件のみ。★**次数 2 以上のものは 1 件も無い。**

```
node tools/absent-recheck.mjs --try 'H2InfRes|infNatTrans_exact|map₂_one|H2_infRes|inflationRestriction'
```
→ `0 件(.cache/mathlib-index.txt、大小無視)`

```
awk -F'\t' '$2 ~ /^groupCohomology\..*(one|One)/ {print $2"\t"$3}' .cache/mathlib-index.txt
```
→ `groupCohomology.map₁_one`(Functoriality.lean:332) は在るが `map₂_one` は無い。

REPL の `#check`:

| 名前 | 結果 |
|---|---|
| `groupCohomology.H1InfRes` / `H1InfRes_exact` | 在る(★次数 1 だけ) |
| `groupCohomology.infNatTrans` | 在る(★全次数。射は引ける) |
| `groupCohomology.H2InfRes` | ★`Unknown constant` |
| `groupCohomology.map₂_one` / `map_one` | ★`Unknown constant` |
| `groupCohomology.mem_coboundaries₂_iff` | ★`Unknown constant`(`coboundaries₂` は `range (d₁₂ A)`) |

★**結論**: `H²` の inflation-restriction は mathlib に無い。本ファイルが埋めた。

## ★★設計 —— 抽象核と具体層

★★**分岐・付値・Galois・コホモロジーの語彙が 1 語も出ない核**を 2 つ切り出した。
どちらも `[Group G] [AddCommGroup A]` と `ρ : G → A →+ A` だけで書かれており、
`k` も `Rep` も `groupCohomology` も出てこない。

* **抽象核 I**(`dOne` / `dOne_sub_isCocycleOn` / `exists_invariant_cochain`) ——
  単射性の中身。「`inf x` が余境界なら、その原始関数 `y` は
  剰余類にしか依らずかつ `S`-不変な `y'` に取り替えられる」。
* **抽象核 II**(`IsCocycle₂` 以下 11 個) —— 完全性の中身。
  「2-コサイクルが `S × S` 上で消えるなら、余境界を引いて
  `S` 方向に完全に潰せる」。★`H¹(S,·) = 0` も
  `hH1 : ∀ f, (1-コサイクル条件) → ∃ a, f s = ρ s a - a` という**語彙ゼロの述語**で受ける。

具体層は `repAddHom`(`Rep k G` の作用を `A →+ A` の族として見る口)を通して核に代入するだけである。

## 逸脱の記録

* **逸脱 1**: 原典 (pGC) は `Γ_K` の**連続**コホモロジーを使うが、本ファイルは
  mathlib の**離散**群コホモロジー `groupCohomology` で書いている。
  ★これは `GroupCohomologyFinite.lean` の逸脱 1 と同じ理由・同じ範囲であり、
  ★本ファイルは「有限商 `G ⧸ S` で得た情報を `G` に上げる」段の**離散版**を与える。
  連続版に上げるには (c) の 2 つ目(有向系の colimit)が要る。
* **逸脱 2**: 原典は 5 項完全列を名指ししていない(「well-known」で畳んでいる)。
  本ファイルは `H¹(S,A) = 0` を仮定に置いた形で右端 2 項の完全性だけを与えている。
  ★transgression 写像 `H¹(S,A)^{G/S} → H²(G/S,A^S)` そのものは**作っていない**。
  仮定を外して 5 項全体を作るには transgression の構成が要る(新ノード)。

## 退化の自己検査

* `inflation₂_bijective_bot` —— `S = ⊥` のとき仮定 `H¹(⊥,A) = 0` は自動で成り立ち、
  さらに `H²(⊥,A) = 0` なので `inf` は**全単射**になる。
  ★定理が空虚でない(仮定を満たす `S` が実在し、しかも結論が `0 = 0` に潰れない)ことの確認。
* `map₂_one` は `map₁_one` の `n = 2` 版であり、`H2InfRes` の `zero` フィールドが
  そこから出る。★`zero` が自明に成り立つ空の短複体ではない。
-/

namespace ABC3.Found.PGC

open CategoryTheory groupCohomology

universe u v

/-! ## 1. ★★抽象核 I —— 1-コチェインの余境界(語彙ゼロ)

★この節には表現論もコホモロジーも分岐も出てこない。
群 `G`、可換群 `A`、`ρ : G → A →+ A` だけである。 -/

section AbstractCoreDescent

variable {G : Type u} {A : Type v} [Group G] [AddCommGroup A]

/-- **1-コチェイン `y : G → A` の余境界** `(dy)(g,h) = ρ(g) y(h) - y(gh) + y(g)`。

★`k` も `Rep` も要らない。`ρ` は加法的自己準同型の族でよい。 -/
def dOne (ρ : G → A →+ A) (y : G → A) (g h : G) : A := ρ g (y h) - y (g * h) + y g

/-- **★抽象核 I-1** —— `dy` が「剰余類にしか依らない」なら、
`s ↦ y(s) - y(1)` は `S` 上の 1-コサイクルである。

★ここで使うのは `hfac`(`dy` が `S` の右移動で不変)と
`hc`(`y 1` が `S`-不変)の 2 つだけ。`S` の正規性すら要らない。 -/
theorem dOne_sub_isCocycleOn (ρ : G → A →+ A) (hone : ∀ a, ρ 1 a = a)
    (S : Subgroup G) (y : G → A)
    (hfac : ∀ g h s t : G, s ∈ S → t ∈ S → dOne ρ y (g * s) (h * t) = dOne ρ y g h)
    (hc : ∀ s ∈ S, ρ s (y 1) = y 1) :
    ∀ s ∈ S, ∀ t ∈ S, y (s * t) - y 1 = ρ s (y t - y 1) + (y s - y 1) := by
  have h11 : dOne ρ y 1 1 = y 1 := by simp [dOne, hone]
  intro s hs t ht
  have h := hfac 1 1 s t hs ht
  rw [one_mul, one_mul, h11] at h
  have hst : ρ s (y t) - y (s * t) + y s = y 1 := h
  have hcs := hc s hs
  rw [map_sub, hcs]
  linear_combination (norm := abel) -hst

/-- **★★抽象核 I-2** —— 単射性の中身。

`dy` が剰余類にしか依らず、かつ `s ↦ y(s) - y(1)` が `S` 上で余境界
(`∃ a, y s - y 1 = ρ s a - a`)ならば、`y` を**同じ余境界を与える** `y'` に取り替えて、
`y'` が

* `S` の右移動で不変(`y'(gs) = y'(g)`)
* 値が `S`-不変(`ρ s (y' g) = y' g`)

を両方満たすようにできる。★これがそのまま「`G ⧸ S` 上の `A^S`-値 1-コチェイン」になる。

★分岐・付値・Galois・コホモロジーの語彙が 1 語も出ない。 -/
theorem exists_invariant_cochain (ρ : G → A →+ A) (hone : ∀ a, ρ 1 a = a)
    (hmul : ∀ g h : G, ∀ a : A, ρ (g * h) a = ρ g (ρ h a))
    (S : Subgroup G) [S.Normal] (y : G → A)
    (hfac : ∀ g h s t : G, s ∈ S → t ∈ S → dOne ρ y (g * s) (h * t) = dOne ρ y g h)
    (a : A) (ha : ∀ s ∈ S, y s - y 1 = ρ s a - a) :
    ∃ y' : G → A, (∀ g h : G, dOne ρ y' g h = dOne ρ y g h) ∧
      (∀ g : G, ∀ s ∈ S, y' (g * s) = y' g) ∧
      (∀ g : G, ∀ s ∈ S, ρ s (y' g) = y' g) := by
  set y' : G → A := fun g => y g - (ρ g a - a) with hy'
  have hd : ∀ g h : G, dOne ρ y' g h = dOne ρ y g h := by
    intro g h; simp only [hy', dOne, map_sub, hmul g h a]; abel
  have hys : ∀ s ∈ S, y' s = y 1 := by
    intro s hs
    have h := ha s hs
    simp only [hy']
    linear_combination (norm := abel) h
  have hy1 : y' 1 = y 1 := hys 1 S.one_mem
  have hcoset : ∀ g : G, ∀ s ∈ S, y' (g * s) = y' g := by
    intro g s hs
    have e1 : dOne ρ y' g s = dOne ρ y' g 1 := by
      rw [hd, hd]
      have := hfac g 1 1 s S.one_mem hs
      rwa [mul_one, one_mul] at this
    simp only [dOne, hys s hs, hy1, mul_one] at e1
    linear_combination (norm := abel) -e1
  refine ⟨y', hd, hcoset, ?_⟩
  intro g s hs
  have hconj : g⁻¹ * s * g ∈ S := Subgroup.Normal.conj_mem' ‹_› s hs g
  have hsg : y' (s * g) = y' g := by
    have h : s * g = g * (g⁻¹ * s * g) := by group
    rw [h, hcoset g _ hconj]
  have e2 : dOne ρ y' s g = dOne ρ y' 1 g := by
    rw [hd, hd]
    have := hfac 1 g s 1 hs S.one_mem
    rwa [one_mul, mul_one] at this
  simp only [dOne, hsg, hys s hs, hy1, one_mul, hone] at e2
  linear_combination (norm := abel) e2

end AbstractCoreDescent

/-! ## 2. ★★抽象核 II —— 2-コサイクル(語彙ゼロ)

★ここも表現論・コホモロジー・分岐の語彙は 1 語も出ない。
「`H¹(S, A) = 0`」ですら `hH1` という**関数についての述語**として受ける。 -/

section AbstractCoreCocycle

variable {G : Type u} {A : Type v} [Group G] [AddCommGroup A]

/-- **2-コサイクル条件**(語彙ゼロ版)。

★mathlib の `groupCohomology.mem_cocycles₂_iff` と**字面まで同じ**形にしてある
(そちらは `G × G → A`、こちらはカリー化)。 -/
def IsCocycle₂ (ρ : G → A →+ A) (x : G → G → A) : Prop :=
  ∀ g h j : G, x (g * h) j + x g h = ρ g (x h j) + x g (h * j)

/-- 1-コチェインの余境界は 2-コサイクル。 -/
theorem isCocycle₂_dOne (ρ : G → A →+ A)
    (hmul : ∀ g h : G, ∀ a : A, ρ (g * h) a = ρ g (ρ h a)) (v : G → A) :
    IsCocycle₂ ρ (dOne ρ v) := by
  intro g h j
  simp only [dOne, map_sub, map_add, hmul g h (v j), mul_assoc]
  abel

/-- 2-コサイクルの差は 2-コサイクル。 -/
theorem IsCocycle₂.sub {ρ : G → A →+ A} {x z : G → G → A}
    (hx : IsCocycle₂ ρ x) (hz : IsCocycle₂ ρ z) :
    IsCocycle₂ ρ (fun g h => x g h - z g h) := by
  intro g h j
  have a1 := hx g h j
  have a2 := hz g h j
  simp only [map_sub]
  linear_combination (norm := abel) a1 - a2

/-- 正規化: `x 1 1 = 0` なら `x 1 g = 0`。 -/
theorem IsCocycle₂.one_left {ρ : G → A →+ A} {x : G → G → A} (hone : ∀ a, ρ 1 a = a)
    (hx : IsCocycle₂ ρ x) (h11 : x 1 1 = 0) (g : G) : x 1 g = 0 := by
  have h := hx 1 1 g
  rw [one_mul, one_mul, h11, hone] at h
  linear_combination (norm := abel) -h

/-- 正規化: `x 1 1 = 0` なら `x g 1 = 0`。 -/
theorem IsCocycle₂.one_right {ρ : G → A →+ A} {x : G → G → A}
    (hx : IsCocycle₂ ρ x) (h11 : x 1 1 = 0) (g : G) : x g 1 = 0 := by
  have h := hx g 1 1
  rw [mul_one, mul_one, h11, map_zero] at h
  linear_combination (norm := abel) h

/-- **★条件 (A) から従う「第 2 引数は剰余類にしか依らない」**。

条件 (A) は `∀ g, ∀ s ∈ S, x g s = 0`。★正規性は要らない。 -/
theorem cocycle₂_right_coset_invariant {ρ : G → A →+ A} {x : G → G → A}
    (hx : IsCocycle₂ ρ x) {S : Subgroup G} (hA : ∀ g : G, ∀ s ∈ S, x g s = 0)
    (g h : G) (s : G) (hs : s ∈ S) : x g (h * s) = x g h := by
  have e := hx g h s
  rw [hA (g * h) s hs, hA h s hs, map_zero] at e
  linear_combination (norm := abel) -e

/-- **★条件 (A)+(B) から「第 1 引数も剰余類にしか依らない」**。

条件 (B) は `∀ s ∈ S, ∀ g, x s g = 0`。★ここで `S` の正規性を使う。 -/
theorem cocycle₂_left_coset_invariant {ρ : G → A →+ A} {x : G → G → A}
    (hx : IsCocycle₂ ρ x) {S : Subgroup G} [S.Normal]
    (hA : ∀ g : G, ∀ s ∈ S, x g s = 0) (hB : ∀ s ∈ S, ∀ g : G, x s g = 0)
    (g h : G) (s : G) (hs : s ∈ S) : x (g * s) h = x g h := by
  have e := hx g s h
  rw [hA g s hs, hB s hs h, map_zero] at e
  have e2 : x (g * s) h = x g (s * h) := by linear_combination (norm := abel) e
  have hconj : h⁻¹ * s * h ∈ S := Subgroup.Normal.conj_mem' ‹_› s hs h
  have e3 : s * h = h * (h⁻¹ * s * h) := by group
  rw [e2, e3, cocycle₂_right_coset_invariant hx hA g h _ hconj]

/-- **★条件 (A)+(B) から値が `S`-不変**。 -/
theorem cocycle₂_invariant_values {ρ : G → A →+ A} {x : G → G → A}
    (hx : IsCocycle₂ ρ x) {S : Subgroup G} [S.Normal]
    (hA : ∀ g : G, ∀ s ∈ S, x g s = 0) (hB : ∀ s ∈ S, ∀ g : G, x s g = 0)
    (g h : G) (u : G) (hu : u ∈ S) : ρ u (x g h) = x g h := by
  have e := hx u g h
  rw [hB u hu g, hB u hu (g * h)] at e
  have hconj : g⁻¹ * u * g ∈ S := Subgroup.Normal.conj_mem' ‹_› u hu g
  have e3 : u * g = g * (g⁻¹ * u * g) := by group
  rw [e3, cocycle₂_left_coset_invariant hx hA hB g h _ hconj] at e
  linear_combination (norm := abel) -e

/-- **★★抽象核 II-1** —— 条件 (A) を余境界で達成する。

`x` が `S × S` 上で消えているなら、代表元関数 `r` から 1-コチェイン `v` を作って
`x - dv` を「第 2 引数が `S` のとき恒等的に 0」にできる。

★`r` は「`(r g)⁻¹ g ∈ S` かつ `r (g s) = r g`」を満たす任意の関数でよい
(すなわち左剰余類の代表元)。★選択公理はここでは使わない(呼ぶ側が `r` を渡す)。 -/
theorem exists_cochain_right {ρ : G → A →+ A}
    {x : G → G → A} (hx : IsCocycle₂ ρ x) {S : Subgroup G}
    (hS : ∀ s ∈ S, ∀ t ∈ S, x s t = 0)
    (r : G → G) (hr1 : ∀ g : G, (r g)⁻¹ * g ∈ S)
    (hr2 : ∀ g : G, ∀ s ∈ S, r (g * s) = r g) :
    ∃ v : G → A, (∀ s ∈ S, v s = 0) ∧ (∀ g : G, ∀ s ∈ S, dOne ρ v g s = x g s) := by
  have hrS : ∀ s ∈ S, r s ∈ S := by
    intro s hs
    have h1 : (r s)⁻¹ * s * s⁻¹ ∈ S := S.mul_mem (hr1 s) (S.inv_mem hs)
    have h2 : (r s)⁻¹ ∈ S := by simpa [mul_assoc] using h1
    simpa using S.inv_mem h2
  have hv0 : ∀ s ∈ S, -x (r s) ((r s)⁻¹ * s) = 0 := by
    intro s hs
    simp [hS _ (hrS s hs) _ (hr1 s)]
  refine ⟨fun g => -x (r g) ((r g)⁻¹ * g), hv0, ?_⟩
  intro g s hs
  have key := hx (r g) ((r g)⁻¹ * g) s
  rw [hS _ (hr1 g) _ hs, map_zero, mul_inv_cancel_left] at key
  show ρ g (-x (r s) ((r s)⁻¹ * s)) -
    (-x (r (g * s)) ((r (g * s))⁻¹ * (g * s))) + (-x (r g) ((r g)⁻¹ * g)) = x g s
  rw [hv0 s hs, map_zero, hr2 g s hs, ← mul_assoc]
  linear_combination (norm := abel) -key

/-- **★★抽象核 II-2** —— 条件 (A) のもとで条件 (B) を余境界で達成する。
★ここが「`H¹(S, A) = 0`」を使う唯一の場所である。

`hH1` は `H¹(S, A) = 0` の**語彙ゼロの言い換え**:
「`S` 上で 1-コサイクル条件を満たす `f : G → A` は `S` 上で余境界である」。

★得られる `v` は (A) を壊さない(`v` が `S` 上で 0 かつ剰余類にしか依らないから)。 -/
theorem exists_cochain_left {ρ : G → A →+ A}
    {x : G → G → A} (hx : IsCocycle₂ ρ x) {S : Subgroup G} [S.Normal]
    (hA : ∀ g : G, ∀ s ∈ S, x g s = 0)
    (hH1 : ∀ f : G → A, (∀ s ∈ S, ∀ t ∈ S, f (s * t) = ρ s (f t) + f s) →
      ∃ a : A, ∀ s ∈ S, f s = ρ s a - a)
    (r : G → G) (hr1 : ∀ g : G, (r g)⁻¹ * g ∈ S)
    (hr2 : ∀ g : G, ∀ s ∈ S, r (g * s) = r g) :
    ∃ v : G → A, (∀ s ∈ S, v s = 0) ∧ (∀ g : G, ∀ s ∈ S, v (g * s) = v g) ∧
      (∀ s ∈ S, ∀ g : G, dOne ρ v s g = x s g) := by
  classical
  have hcocL : ∀ g : G, ∀ s ∈ S, ∀ t ∈ S, x (s * t) g = ρ s (x t g) + x s g := by
    intro g s hs t ht
    have e := hx s t g
    rw [hA s t ht] at e
    have hconj : g⁻¹ * t * g ∈ S := Subgroup.Normal.conj_mem' ‹_› t ht g
    have e3 : t * g = g * (g⁻¹ * t * g) := by group
    rw [e3, cocycle₂_right_coset_invariant hx hA s g _ hconj] at e
    linear_combination (norm := abel) e
  choose a ha using fun g => hH1 (fun s => x s g) (hcocL g)
  refine ⟨fun g => if g ∈ S then 0 else a (r g), ?_, ?_, ?_⟩
  · intro s hs; simp [hs]
  · intro g s hs
    by_cases hg : g ∈ S
    · simp [hg, S.mul_mem hg hs]
    · have hgs : g * s ∉ S := by
        intro hc
        exact hg (by simpa [mul_assoc] using S.mul_mem hc (S.inv_mem hs))
      simp [hg, hgs, hr2 g s hs]
  · intro s hs g
    have hconj : g⁻¹ * s * g ∈ S := Subgroup.Normal.conj_mem' ‹_› s hs g
    have e3 : s * g = g * (g⁻¹ * s * g) := by group
    have hvsg : (if s * g ∈ S then (0 : A) else a (r (s * g)))
        = (if g ∈ S then (0 : A) else a (r g)) := by
      rw [e3]
      by_cases hg : g ∈ S
      · simp [hg, S.mul_mem hg hconj]
      · have hgs : g * (g⁻¹ * s * g) ∉ S := by
          intro hc
          exact hg (by simpa [mul_assoc] using S.mul_mem hc (S.inv_mem hconj))
        simp [hg, hgs, hr2 g _ hconj]
    show ρ s (if g ∈ S then (0 : A) else a (r g)) - (if s * g ∈ S then (0 : A) else a (r (s * g)))
      + (if s ∈ S then (0 : A) else a (r s)) = x s g
    rw [hvsg, if_pos hs, add_zero]
    by_cases hg : g ∈ S
    · rw [if_pos hg, map_zero, sub_zero, hA s g hg]
    · rw [if_neg hg]
      have h1 : x s (r g) = ρ s (a (r g)) - a (r g) := ha (r g) s hs
      have h2 : x s g = x s (r g) := by
        conv_lhs => rw [show g = r g * ((r g)⁻¹ * g) from (mul_inv_cancel_left _ _).symm]
        exact cocycle₂_right_coset_invariant hx hA s (r g) _ (hr1 g)
      rw [h2, h1]

/-- **★★★抽象核 II-3(完全性の中身)** —— `S × S` 上で消える 2-コサイクルは、
余境界を引くと「剰余類にしか依らず、値が `S`-不変」になる。

★すなわち `x - dv` は `G ⧸ S` 上の `A^S`-値 2-コサイクルから来ている。
これが inflation-restriction の完全性の核である。

★分岐・付値・Galois・コホモロジーの語彙が 1 語も出ない。 -/
theorem exists_inflation_cochain {ρ : G → A →+ A}
    (hmul : ∀ g h : G, ∀ a : A, ρ (g * h) a = ρ g (ρ h a))
    {x : G → G → A} (hx : IsCocycle₂ ρ x) {S : Subgroup G} [S.Normal]
    (hS : ∀ s ∈ S, ∀ t ∈ S, x s t = 0)
    (hH1 : ∀ f : G → A, (∀ s ∈ S, ∀ t ∈ S, f (s * t) = ρ s (f t) + f s) →
      ∃ a : A, ∀ s ∈ S, f s = ρ s a - a)
    (r : G → G) (hr1 : ∀ g : G, (r g)⁻¹ * g ∈ S)
    (hr2 : ∀ g : G, ∀ s ∈ S, r (g * s) = r g) :
    ∃ v : G → A,
      (∀ g h : G, ∀ s ∈ S, ∀ t ∈ S,
        x (g * s) (h * t) - dOne ρ v (g * s) (h * t) = x g h - dOne ρ v g h) ∧
      (∀ g h : G, ∀ u ∈ S, ρ u (x g h - dOne ρ v g h) = x g h - dOne ρ v g h) := by
  obtain ⟨v₀, hv₀0, hv₀⟩ := exists_cochain_right hx hS r hr1 hr2
  have hx₁ : IsCocycle₂ ρ (fun g h => x g h - dOne ρ v₀ g h) :=
    hx.sub (isCocycle₂_dOne ρ hmul v₀)
  have hA₁ : ∀ g : G, ∀ s ∈ S, (fun g h => x g h - dOne ρ v₀ g h) g s = 0 := by
    intro g s hs
    show x g s - dOne ρ v₀ g s = 0
    rw [hv₀ g s hs]; abel
  obtain ⟨v₁, hv₁0, hv₁c, hv₁⟩ := exists_cochain_left hx₁ hA₁ hH1 r hr1 hr2
  have hsplit : ∀ g h : G, x g h - dOne ρ (fun g => v₀ g + v₁ g) g h
      = (fun g h => x g h - dOne ρ v₀ g h) g h - dOne ρ v₁ g h := by
    intro g h; simp only [dOne, map_add]; abel
  have hA₂ : ∀ g : G, ∀ s ∈ S,
      (fun g h => (fun g h => x g h - dOne ρ v₀ g h) g h - dOne ρ v₁ g h) g s = 0 := by
    intro g s hs
    have d1 : dOne ρ v₁ g s = 0 := by
      simp only [dOne, hv₁0 s hs, hv₁c g s hs, map_zero]; abel
    show x g s - dOne ρ v₀ g s - dOne ρ v₁ g s = 0
    rw [hv₀ g s hs, d1]; abel
  have hB₂ : ∀ s ∈ S, ∀ g : G,
      (fun g h => (fun g h => x g h - dOne ρ v₀ g h) g h - dOne ρ v₁ g h) s g = 0 := by
    intro s hs g
    show x s g - dOne ρ v₀ s g - dOne ρ v₁ s g = 0
    rw [hv₁ s hs g]; abel
  have hx₂ : IsCocycle₂ ρ (fun g h => (fun g h => x g h - dOne ρ v₀ g h) g h - dOne ρ v₁ g h) :=
    hx₁.sub (isCocycle₂_dOne ρ hmul v₁)
  refine ⟨fun g => v₀ g + v₁ g, ?_, ?_⟩
  · intro g h s hs t ht
    rw [hsplit, hsplit]
    rw [cocycle₂_right_coset_invariant hx₂ hA₂ _ _ _ ht,
      cocycle₂_left_coset_invariant hx₂ hA₂ hB₂ _ _ _ hs]
  · intro g h u hu
    rw [hsplit]
    exact cocycle₂_invariant_values hx₂ hA₂ hB₂ g h u hu

end AbstractCoreCocycle

/-! ## 3. 具体層への口 —— `Rep k G` と抽象核の橋

★ここから `Rep` が出てくる。核に代入するための最小限の配管だけ。 -/

section Bridge

variable {k G : Type} [CommRing k] [Group G]

/-- **`Rep k G` の作用を加法的自己準同型の族として見る口**(抽象核へ渡すため)。 -/
noncomputable def repAddHom (A : Rep k G) (g : G) : A →+ A := (A.ρ g).toAddMonoidHom

theorem repAddHom_apply (A : Rep k G) (g : G) (a : A) : repAddHom A g a = A.ρ g a := rfl

theorem repAddHom_one (A : Rep k G) (a : A) : repAddHom A 1 a = a := by
  simp [repAddHom_apply]

theorem repAddHom_mul (A : Rep k G) (g h : G) (a : A) :
    repAddHom A (g * h) a = repAddHom A g (repAddHom A h a) := by
  simp [repAddHom_apply]

/-- **左剰余類の代表元**(選択公理)。抽象核 II の `r` に渡すためのもの。 -/
noncomputable def cosetRep (S : Subgroup G) (g : G) : G := (QuotientGroup.mk g : G ⧸ S).out

theorem cosetRep_inv_mul_mem (S : Subgroup G) (g : G) : (cosetRep S g)⁻¹ * g ∈ S :=
  QuotientGroup.eq.1 (Quotient.out_eq _)

theorem cosetRep_mul_of_mem (S : Subgroup G) (g s : G) (hs : s ∈ S) :
    cosetRep S (g * s) = cosetRep S g := by
  unfold cosetRep
  rw [QuotientGroup.mk_mul_of_mem g hs]

/-- **`H¹(S, A) = 0` を抽象核が要求する形に直す**。

★`Subsingleton ↑(groupCohomology (Rep.res S.subtype A) 1)` が本ファイル唯一の仮定である。 -/
theorem hH1_of_subsingleton (A : Rep k G) (S : Subgroup G)
    (h1 : Subsingleton (groupCohomology (Rep.res S.subtype A) 1))
    (f : G → A) (hf : ∀ s ∈ S, ∀ t ∈ S, f (s * t) = repAddHom A s (f t) + f s) :
    ∃ a : A, ∀ s ∈ S, f s = repAddHom A s a - a := by
  have hmem : (fun s : S => f (s : G)) ∈ cocycles₁ (Rep.res S.subtype A) := by
    rw [mem_cocycles₁_iff]
    intro g h
    exact hf (g : G) g.2 (h : G) h.2
  have h0 : (ConcreteCategory.hom (H1π (Rep.res S.subtype A)))
      ⟨(fun s : S => f (s : G)), hmem⟩ = 0 := h1.elim _ _
  obtain ⟨a, ha⟩ := (H1π_eq_zero_iff _).1 h0
  exact ⟨a, fun s hs => (congrFun ha ⟨s, hs⟩).symm⟩

end Bridge

/-! ## 4. ★★mathlib の欠落 `map₂_one` と短複体 `H2InfRes`

mathlib は `groupCohomology.map₁_one`(`map 1 φ 1 = 0`)しか持たない。
`H2InfRes` の `zero` フィールドにはその `n = 2` 版が要る。 -/

section ShortComplexDef

variable {k G : Type} [CommRing k] [Group G]

/-- **★mathlib に無い**: 自明な準同型に沿った `H²` の写像は 0。

★`map₁_one`(Functoriality.lean:332)の `n = 2` 版。
証明の要点は「`φ : Res(1)(A) ⟶ B` の像は必ず `G`-不変」であること。
そこから `mapCocycles₂ 1 φ x` が**定数** `φ(x(1,1))` になり、
定数 `S`-不変コチェインは自分自身の余境界なので `H²` では 0 になる。 -/
theorem map₂_one {H : Type} [Group H] {A : Rep k H} {B : Rep k G}
    (φ : Rep.res (1 : G →* H) A ⟶ B) :
    groupCohomology.map (1 : G →* H) φ 2 = 0 := by
  refine ModuleCat.hom_ext (LinearMap.ext fun z => ?_)
  induction z using groupCohomology.H2_induction_on with
  | h x =>
    show (ConcreteCategory.hom (groupCohomology.map (1 : G →* H) φ 2))
      ((ConcreteCategory.hom (H2π A)) x) = 0
    rw [H2π_comp_map_apply, H2π_eq_zero_iff]
    refine ⟨fun _ => φ.hom (x (1, 1)), ?_⟩
    funext p
    have hinv : ∀ g : G, B.ρ g (φ.hom (x (1, 1))) = φ.hom (x (1, 1)) := by
      intro g
      have := Rep.hom_comm_apply φ g (x (1, 1))
      simpa using this.symm
    show B.ρ p.1 (φ.hom (x (1, 1))) - φ.hom (x (1, 1)) + φ.hom (x (1, 1)) = _
    rw [hinv, sub_add_cancel]
    exact (congrFun (coe_mapCocycles₂ (1 : G →* H) φ x) p).symm

/-- **★★短複体 `H²(G ⧸ S, Aˢ) ⟶ H²(G, A) ⟶ H²(S, A)`**。

★mathlib の `groupCohomology.H1InfRes` の `n = 2` 版。
`zero` フィールドは `map₂_one` から出る。 -/
@[simps X₁ X₂ X₃ f g]
noncomputable def H2InfRes (A : Rep k G) (S : Subgroup G) [S.Normal] :
    ShortComplex (ModuleCat k) where
  X₁ := groupCohomology (A.quotientToInvariants S) 2
  X₂ := groupCohomology A 2
  X₃ := groupCohomology (Rep.res S.subtype A) 2
  f := map (QuotientGroup.mk' S) (Rep.ofHom <| A.ρ.quotientToInvariants_lift S) 2
  g := map S.subtype (𝟙 _) 2
  zero := by
    rw [← map_comp, Category.comp_id,
      congr (QuotientGroup.mk'_comp_subtype S) (fun f φ => map f φ 2), map₂_one]

end ShortComplexDef

/-! ## 5. ★★★inflation の単射性(= 5 項完全列の右端の完全性)

★古典的な 5 項完全列
`H¹(G/S,A^S) → H¹(G,A) → H¹(S,A)^{G/S} → H²(G/S,A^S) → H²(G,A)`
の右端 `H²(G/S,A^S)` における完全性は、`H¹(S,A) = 0` のとき
「inflation が単射」に一致する(transgression の source が消えるから)。 -/

section Injectivity

variable {k G : Type} [CommRing k] [Group G]

def inflation₂_injective.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 3, item := "Proposition 1.1", sectionId := "prop-1-1" }

/-- **★コサイクルの水準での単射性**(仮定 `H¹(S,A) = 0` を使わない形)。

`x` を `G ⧸ S` 上の 2-コサイクルとし、その inflation が `y : G → A` の余境界であるとする。
もし `s ↦ y s - y 1` が `S` 上で余境界(`∃ a`)なら、`x` 自身が `G ⧸ S` 上で余境界である。

★中身はすべて抽象核 I-2(`exists_invariant_cochain`)。ここは `Rep` への翻訳だけ。 -/
theorem H2π_eq_zero_of_split (A : Rep k G) (S : Subgroup G) [S.Normal]
    (x : cocycles₂ (A.quotientToInvariants S)) (y : G → A)
    (hy : ∀ g h : G, dOne (repAddHom A) y g h
        = (x (QuotientGroup.mk' S g, QuotientGroup.mk' S h) : A))
    (a : A) (ha : ∀ s ∈ S, y s - y 1 = A.ρ s a - a) :
    (ConcreteCategory.hom (H2π (A.quotientToInvariants S))) x = 0 := by
  have hfac : ∀ g h s t : G, s ∈ S → t ∈ S →
      dOne (repAddHom A) y (g * s) (h * t) = dOne (repAddHom A) y g h := by
    intro g h s t hs ht
    rw [hy, hy]
    have e1 : QuotientGroup.mk' S (g * s) = QuotientGroup.mk' S g := by
      simpa using QuotientGroup.mk_mul_of_mem (s := S) g hs
    have e2 : QuotientGroup.mk' S (h * t) = QuotientGroup.mk' S h := by
      simpa using QuotientGroup.mk_mul_of_mem (s := S) h ht
    rw [e1, e2]
  obtain ⟨y', hd, hcoset, hinvS⟩ :=
    exists_invariant_cochain (repAddHom A) (repAddHom_one A) (repAddHom_mul A) S y hfac a ha
  have hmem : ∀ g : G, y' g ∈ Representation.invariants (A.ρ.comp S.subtype) := fun g =>
    (Representation.mem_invariants _ _).2 (fun s => hinvS g s s.2)
  have hwd : ∀ g h : G, g⁻¹ * h ∈ S → y' g = y' h := by
    intro g h hs
    have e : h = g * (g⁻¹ * h) := by group
    rw [e, hcoset g _ hs]
  rw [H2π_eq_zero_iff]
  refine ⟨fun q => Quotient.liftOn' q (fun g => (⟨y' g, hmem g⟩ : (A.quotientToInvariants S)))
    (fun g h hgh => Subtype.ext (hwd g h (QuotientGroup.leftRel_apply.1 hgh))), ?_⟩
  funext p
  obtain ⟨q1, q2⟩ := p
  induction q1 using QuotientGroup.induction_on with | @H g =>
  induction q2 using QuotientGroup.induction_on with | @H h =>
  apply Subtype.ext
  show dOne (repAddHom A) y' g h = (x (QuotientGroup.mk' S g, QuotientGroup.mk' S h) : A)
  exact (hd g h).trans (hy g h)

/-- **★★★`H¹(S, A) = 0` ならば inflation `H²(G ⧸ S, Aˢ) ⟶ H²(G, A)` は単射**。

原文 (pGC p.3):
> This is clearly a group-theoretic condition on M.

★mathlib は `H1InfRes` の `Mono` インスタンスを持つが、`n = 2` は持たない
(`absent-recheck --try 'H2InfRes|...'` → 0 件)。

★仮定は `H¹(S, A) = 0` のみ。★これは 5 項完全列の右端の完全性そのものである。 -/
theorem inflation₂_injective (A : Rep k G) (S : Subgroup G) [S.Normal]
    (h1 : Subsingleton (groupCohomology (Rep.res S.subtype A) 1)) :
    Function.Injective (ConcreteCategory.hom
      (map (A := A.quotientToInvariants S) (B := A) (QuotientGroup.mk' S)
        (Rep.ofHom (A.ρ.quotientToInvariants_lift S)) 2)) := by
  rw [injective_iff_map_eq_zero]
  intro z hz
  induction z using groupCohomology.H2_induction_on with
  | @h x =>
  rw [H2π_comp_map_apply, H2π_eq_zero_iff] at hz
  obtain ⟨y, hy0⟩ := hz
  have hy : ∀ g h : G, dOne (repAddHom A) y g h
      = (x (QuotientGroup.mk' S g, QuotientGroup.mk' S h) : A) := fun g h => congrFun hy0 (g, h)
  have h11 : dOne (repAddHom A) y 1 1 = y 1 := by simp [dOne, repAddHom_one]
  have hx11 : y 1 = ((x (1, 1) : (A.quotientToInvariants S)) : A) := by
    have h := h11.symm.trans (hy 1 1)
    simpa using h
  have hc : ∀ s ∈ S, repAddHom A s (y 1) = y 1 := by
    intro s hs
    rw [hx11, repAddHom_apply]
    exact (Representation.mem_invariants _ _).1 (x (1, 1)).2 ⟨s, hs⟩
  have hfac : ∀ g h s t : G, s ∈ S → t ∈ S →
      dOne (repAddHom A) y (g * s) (h * t) = dOne (repAddHom A) y g h := by
    intro g h s t hs ht
    rw [hy, hy]
    have e1 : QuotientGroup.mk' S (g * s) = QuotientGroup.mk' S g := by
      simpa using QuotientGroup.mk_mul_of_mem (s := S) g hs
    have e2 : QuotientGroup.mk' S (h * t) = QuotientGroup.mk' S h := by
      simpa using QuotientGroup.mk_mul_of_mem (s := S) h ht
    rw [e1, e2]
  have hcoc := dOne_sub_isCocycleOn (repAddHom A) (repAddHom_one A) S y hfac hc
  have humem : (fun s : S => y (s : G) - y 1) ∈ cocycles₁ (Rep.res S.subtype A) := by
    rw [mem_cocycles₁_iff]
    intro g h
    exact hcoc (g : G) g.2 (h : G) h.2
  have hu0 : (ConcreteCategory.hom (H1π (Rep.res S.subtype A)))
      ⟨(fun s : S => y (s : G) - y 1), humem⟩ = 0 := h1.elim _ _
  obtain ⟨a, ha0⟩ := (H1π_eq_zero_iff _).1 hu0
  have ha : ∀ s ∈ S, y s - y 1 = A.ρ s a - a := by
    intro s hs
    exact (congrFun ha0 ⟨s, hs⟩).symm
  exact H2π_eq_zero_of_split A S x y hy a ha

/-- 短複体 `H2InfRes` の左端が単射。 -/
theorem mono_H2InfRes_f (A : Rep k G) (S : Subgroup G) [S.Normal]
    (h1 : Subsingleton (groupCohomology (Rep.res S.subtype A) 1)) :
    Mono (H2InfRes A S).f :=
  (ModuleCat.mono_iff_injective _).2 (inflation₂_injective A S h1)

/-- **★直前の波が作った `inflationH2` がそのまま単射になる**
(`GroupCohomologyFinite.lean` の `inflationH2` = `infNatTrans` の `n = 2` 評価)。 -/
theorem inflationH2_injective (A : Rep k G) (S : Subgroup G) [S.Normal]
    (h1 : Subsingleton (groupCohomology (Rep.res S.subtype A) 1)) :
    Function.Injective (ConcreteCategory.hom (inflationH2 A S)) :=
  inflation₂_injective A S h1

theorem mono_inflationH2 (A : Rep k G) (S : Subgroup G) [S.Normal]
    (h1 : Subsingleton (groupCohomology (Rep.res S.subtype A) 1)) :
    Mono (inflationH2 A S) :=
  (ModuleCat.mono_iff_injective _).2 (inflationH2_injective A S h1)

/-- `IsZero` 版(圏論的な言い方)。 -/
theorem inflationH2_injective_of_isZero (A : Rep k G) (S : Subgroup G) [S.Normal]
    (h1 : Limits.IsZero (groupCohomology (Rep.res S.subtype A) 1)) :
    Function.Injective (ConcreteCategory.hom (inflationH2 A S)) :=
  inflationH2_injective A S (ModuleCat.isZero_iff_subsingleton.mp h1)

end Injectivity

/-! ## 6. ★★★`H²(G, A)` における完全性

`H¹(S,A) = 0` のもとで `H²(G/S, A^S) → H²(G,A) → H²(S,A)` は完全。
★古典的な「`H^i(S,A) = 0` for `i < n` ならば
`0 → H^n(G/S,A^S) → H^n(G,A) → H^n(S,A)` が完全」の `n = 2` の場合。 -/

section Exactness

variable {k G : Type} [CommRing k] [Group G]

/-- 抽象核 II-3 が出した「剰余類にしか依らず `S`-不変な」関数を
`G ⧸ S` 上の `Aˢ`-値 2-コチェインに降ろす。 -/
noncomputable def descendCocycle (A : Rep k G) (S : Subgroup G) [S.Normal] (x₂ : G → G → A)
    (hmem : ∀ g h : G, x₂ g h ∈ Representation.invariants (A.ρ.comp S.subtype))
    (hfac : ∀ g h : G, ∀ s ∈ S, ∀ t ∈ S, x₂ (g * s) (h * t) = x₂ g h) :
    (G ⧸ S) × (G ⧸ S) → (A.quotientToInvariants S) := fun p =>
  Quotient.liftOn₂' p.1 p.2 (fun g h => (⟨x₂ g h, hmem g h⟩ : (A.quotientToInvariants S)))
    (by
      intro g1 h1 g2 h2 e1 e2
      have m1 : g1⁻¹ * g2 ∈ S := QuotientGroup.leftRel_apply.1 e1
      have m2 : h1⁻¹ * h2 ∈ S := QuotientGroup.leftRel_apply.1 e2
      apply Subtype.ext
      show x₂ g1 h1 = x₂ g2 h2
      conv_rhs => rw [show g2 = g1 * (g1⁻¹ * g2) from (mul_inv_cancel_left _ _).symm,
        show h2 = h1 * (h1⁻¹ * h2) from (mul_inv_cancel_left _ _).symm]
      exact (hfac g1 h1 _ m1 _ m2).symm)

theorem descendCocycle_mk (A : Rep k G) (S : Subgroup G) [S.Normal] (x₂ : G → G → A)
    (hmem : ∀ g h : G, x₂ g h ∈ Representation.invariants (A.ρ.comp S.subtype))
    (hfac : ∀ g h : G, ∀ s ∈ S, ∀ t ∈ S, x₂ (g * s) (h * t) = x₂ g h) (g h : G) :
    ((descendCocycle A S x₂ hmem hfac (QuotientGroup.mk' S g, QuotientGroup.mk' S h) : _) : A)
      = x₂ g h := rfl

theorem descendCocycle_mem (A : Rep k G) (S : Subgroup G) [S.Normal] (x₂ : G → G → A)
    (hmem : ∀ g h : G, x₂ g h ∈ Representation.invariants (A.ρ.comp S.subtype))
    (hfac : ∀ g h : G, ∀ s ∈ S, ∀ t ∈ S, x₂ (g * s) (h * t) = x₂ g h)
    (hcoc : IsCocycle₂ (repAddHom A) x₂) :
    descendCocycle A S x₂ hmem hfac ∈ cocycles₂ (A.quotientToInvariants S) := by
  rw [mem_cocycles₂_iff]
  intro q1 q2 q3
  induction q1 using QuotientGroup.induction_on with | @H g =>
  induction q2 using QuotientGroup.induction_on with | @H h =>
  induction q3 using QuotientGroup.induction_on with | @H j =>
  apply Subtype.ext
  show x₂ (g * h) j + x₂ g h = repAddHom A g (x₂ h j) + x₂ g (h * j)
  exact hcoc g h j

def exists_inflation_preimage.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 3, item := "Proposition 1.1", sectionId := "prop-1-1" }

/-- **★★★`H¹(S, A) = 0` のとき、`H²(S, A)` へ 0 に行く類は inflation の像である**。

原文 (pGC p.3):
> Thus, we conclude that the isomorphism class of the ΓK-module Zp(1) can be recovered
> group-theoretically from ΓK.

★段取り(すべて抽象核 II にある):
1. `res` が 0 なので `z` は `S × S` 上で余境界 `dw` に一致する。`w` を `G` へ延長して引く。
2. 抽象核 II-1 で「第 2 引数が `S` なら 0」(条件 (A))にする。
3. 抽象核 II-2(★ここだけが `H¹(S,A) = 0` を使う)で「第 1 引数が `S` なら 0」(条件 (B))にする。
4. 抽象核 II-3 の結論より、残りは剰余類にしか依らず値が `S`-不変。`G ⧸ S` へ降ろす。 -/
theorem exists_inflation_preimage (A : Rep k G) (S : Subgroup G) [S.Normal]
    (h1 : Subsingleton (groupCohomology (Rep.res S.subtype A) 1))
    (c : groupCohomology A 2)
    (hc : (ConcreteCategory.hom (map (A := A) (B := Rep.res S.subtype A) S.subtype
        (𝟙 (Rep.res S.subtype A)) 2)) c = 0) :
    ∃ b : groupCohomology (A.quotientToInvariants S) 2,
      (ConcreteCategory.hom (map (A := A.quotientToInvariants S) (B := A) (QuotientGroup.mk' S)
        (Rep.ofHom (A.ρ.quotientToInvariants_lift S)) 2)) b = c := by
  classical
  induction c using groupCohomology.H2_induction_on with
  | @h z =>
  rw [H2π_comp_map_apply, H2π_eq_zero_iff] at hc
  obtain ⟨w, hw⟩ := hc
  have hw' : ∀ s t : S, repAddHom A (s : G) (w t) - w (s * t) + w s = z ((s : G), (t : G)) :=
    fun s t => congrFun hw (s, t)
  set wt : G → A := fun g => if h : g ∈ S then w ⟨g, h⟩ else 0 with hwt
  have hzcoc : IsCocycle₂ (repAddHom A) (fun g h => z (g, h)) := (mem_cocycles₂_iff _).1 z.2
  have hxcoc : IsCocycle₂ (repAddHom A) (fun g h => z (g, h) - dOne (repAddHom A) wt g h) :=
    hzcoc.sub (isCocycle₂_dOne _ (repAddHom_mul A) wt)
  have hSx : ∀ s ∈ S, ∀ t ∈ S,
      (fun g h => z (g, h) - dOne (repAddHom A) wt g h) s t = 0 := by
    intro s hs t ht
    have h : repAddHom A s (w ⟨t, ht⟩) - w ⟨s * t, S.mul_mem hs ht⟩ + w ⟨s, hs⟩ = z (s, t) :=
      hw' ⟨s, hs⟩ ⟨t, ht⟩
    show z (s, t) - (repAddHom A s (wt t) - wt (s * t) + wt s) = 0
    simp only [hwt, dif_pos hs, dif_pos ht, dif_pos (S.mul_mem hs ht)]
    rw [← h]; abel
  obtain ⟨v, hfac, hinv⟩ := exists_inflation_cochain (repAddHom_mul A) hxcoc hSx
    (hH1_of_subsingleton A S h1) (cosetRep S) (cosetRep_inv_mul_mem S) (cosetRep_mul_of_mem S)
  set V : G → A := fun g => wt g + v g with hV
  have hx₂eq : ∀ g h : G, z (g, h) - dOne (repAddHom A) V g h
      = (fun g h => z (g, h) - dOne (repAddHom A) wt g h) g h - dOne (repAddHom A) v g h := by
    intro g h; simp only [hV, dOne, map_add]; abel
  have hfac₂ : ∀ g h : G, ∀ s ∈ S, ∀ t ∈ S,
      (fun g h => z (g, h) - dOne (repAddHom A) V g h) (g * s) (h * t)
        = (fun g h => z (g, h) - dOne (repAddHom A) V g h) g h := by
    intro g h s hs t ht
    show z (g * s, h * t) - dOne (repAddHom A) V (g * s) (h * t)
      = z (g, h) - dOne (repAddHom A) V g h
    rw [hx₂eq, hx₂eq]
    exact hfac g h s hs t ht
  have hinv₂ : ∀ g h : G, ∀ u ∈ S,
      repAddHom A u ((fun g h => z (g, h) - dOne (repAddHom A) V g h) g h)
        = (fun g h => z (g, h) - dOne (repAddHom A) V g h) g h := by
    intro g h u hu
    show repAddHom A u (z (g, h) - dOne (repAddHom A) V g h)
      = z (g, h) - dOne (repAddHom A) V g h
    rw [hx₂eq]
    exact hinv g h u hu
  have hmem₂ : ∀ g h : G, (fun g h => z (g, h) - dOne (repAddHom A) V g h) g h
      ∈ Representation.invariants (A.ρ.comp S.subtype) := fun g h =>
    (Representation.mem_invariants _ _).2 (fun u => hinv₂ g h u u.2)
  have hcoc₂ : IsCocycle₂ (repAddHom A) (fun g h => z (g, h) - dOne (repAddHom A) V g h) :=
    hzcoc.sub (isCocycle₂_dOne _ (repAddHom_mul A) V)
  refine ⟨(ConcreteCategory.hom (H2π (A.quotientToInvariants S)))
    ⟨descendCocycle A S _ hmem₂ hfac₂, descendCocycle_mem A S _ hmem₂ hfac₂ hcoc₂⟩, ?_⟩
  rw [H2π_comp_map_apply, H2π_eq_iff]
  refine ⟨fun g => -V g, ?_⟩
  funext p
  obtain ⟨g, h⟩ := p
  show repAddHom A g (-V h) - (-V (g * h)) + (-V g)
    = ((descendCocycle A S _ hmem₂ hfac₂ (QuotientGroup.mk' S g, QuotientGroup.mk' S h) : _) : A)
      - z (g, h)
  rw [descendCocycle_mk]
  show repAddHom A g (-V h) - (-V (g * h)) + (-V g)
    = (z (g, h) - dOne (repAddHom A) V g h) - z (g, h)
  simp only [dOne, map_neg]
  abel

/-- **★★★短複体 `H²(G ⧸ S, Aˢ) ⟶ H²(G, A) ⟶ H²(S, A)` は `H²(G,A)` において完全**。

★仮定は `H¹(S, A) = 0` のみ。★`mono_H2InfRes_f` と合わせて
`0 → H²(G ⧸ S, Aˢ) → H²(G, A) → H²(S, A)` が完全であることを意味する。 -/
theorem H2InfRes_exact (A : Rep k G) (S : Subgroup G) [S.Normal]
    (h1 : Subsingleton (groupCohomology (Rep.res S.subtype A) 1)) :
    (H2InfRes A S).Exact := by
  rw [ShortComplex.moduleCat_exact_iff_ker_sub_range]
  intro c hc
  exact exists_inflation_preimage A S h1 c hc

/-- `IsZero` 版。 -/
theorem H2InfRes_exact_of_isZero (A : Rep k G) (S : Subgroup G) [S.Normal]
    (h1 : Limits.IsZero (groupCohomology (Rep.res S.subtype A) 1)) :
    (H2InfRes A S).Exact :=
  H2InfRes_exact A S (ModuleCat.isZero_iff_subsingleton.mp h1)

/-- **★★★本ファイルの到達点をひとまとめにしたもの**。

`H¹(S, A) = 0` ならば `0 → H²(G ⧸ S, Aˢ) → H²(G, A) → H²(S, A)` は完全。 -/
theorem mono_and_exact_H2InfRes (A : Rep k G) (S : Subgroup G) [S.Normal]
    (h1 : Subsingleton (groupCohomology (Rep.res S.subtype A) 1)) :
    Mono (H2InfRes A S).f ∧ (H2InfRes A S).Exact :=
  ⟨mono_H2InfRes_f A S h1, H2InfRes_exact A S h1⟩

end Exactness

/-! ## 7. ★退化の自己検査 —— 仮定は空虚ではない

`S = ⊥` のとき仮定 `H¹(⊥, A) = 0` は**自動で成り立つ**。
しかも `H²(⊥, A) = 0` なので、そのとき inflation は**全単射**になる。
★「仮定を満たす `S` が実在し、結論も `0 = 0` に潰れない」ことの確認。 -/

section Degeneracy

variable {k G : Type} [CommRing k] [Group G]

/-- 自明群の `H¹` は 0。 -/
theorem subsingleton_H1_of_subsingleton (B : Rep k G) [Subsingleton G] :
    Subsingleton (groupCohomology B 1) := by
  have key : ∀ x : cocycles₁ B, (x : G → B) = 0 := by
    intro x
    funext g
    rw [Subsingleton.elim g 1]
    simp
  constructor
  intro a b
  induction a using groupCohomology.H1_induction_on with | @h u =>
  induction b using groupCohomology.H1_induction_on with | @h u' =>
  have e : u = u' := Subtype.ext ((key u).trans (key u').symm)
  rw [e]

/-- 自明群の `H²` は 0。★定数コチェインが自分自身の余境界であることによる。 -/
theorem subsingleton_H2_of_subsingleton (B : Rep k G) [Subsingleton G] :
    Subsingleton (groupCohomology B 2) := by
  constructor
  intro a b
  induction a using groupCohomology.H2_induction_on with | @h u =>
  induction b using groupCohomology.H2_induction_on with | @h u' =>
  rw [H2π_eq_iff]
  refine ⟨fun _ => (u : G × G → B) (1, 1) - (u' : G × G → B) (1, 1), ?_⟩
  funext p
  obtain ⟨g, h⟩ := p
  show B.ρ g ((u : G × G → B) (1, 1) - (u' : G × G → B) (1, 1)) - _ + _ = _
  rw [Subsingleton.elim g 1, Subsingleton.elim h 1]
  simp

/-- **★★退化の自己検査** —— `S = ⊥` のとき inflation は全単射。

★単射性は `inflation₂_injective`(仮定は `H¹(⊥,A) = 0`、これは自動)、
全射性は `exists_inflation_preimage`(`H²(⊥,A) = 0` なので前提が自動)から出る。
★これで「仮定 `H¹(S,A) = 0` を満たす `S` が実在し、しかも結論が空虚でない」ことが見える。 -/
theorem inflation₂_bijective_bot (A : Rep k G) :
    Function.Bijective (ConcreteCategory.hom
      (map (A := A.quotientToInvariants (⊥ : Subgroup G)) (B := A)
        (QuotientGroup.mk' (⊥ : Subgroup G))
        (Rep.ofHom (A.ρ.quotientToInvariants_lift (⊥ : Subgroup G))) 2)) := by
  have h1 : Subsingleton (groupCohomology (Rep.res (⊥ : Subgroup G).subtype A) 1) :=
    subsingleton_H1_of_subsingleton _
  refine ⟨inflation₂_injective A ⊥ h1, fun c => ?_⟩
  have h2 : Subsingleton (groupCohomology (Rep.res (⊥ : Subgroup G).subtype A) 2) :=
    subsingleton_H2_of_subsingleton _
  exact exists_inflation_preimage A ⊥ h1 c (h2.elim _ _)

end Degeneracy

end ABC3.Found.PGC
