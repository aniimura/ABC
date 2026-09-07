import ABC3.Found.PGC.SenJumpExpansion

/-!
# 商で測った `i` は剰余類上の平均(Yoshida 2008 Lemma 6.8)

典拠: T. Yoshida, *Local Class Field Theory via Lubin-Tate Theory*
(Ann. Fac. Sci. Toulouse 17-2, 2008; arXiv math/0606108) **Lemma 6.8**(物理 p.15)。
構造化済み原文は `ResearchPaper/1_Structured/Local Class Field Theory via Lubin-Tate Theory/`
`section-6.html` の `id="lemma-6-8"`(`data-pdf-page="15"`, `data-item="Lemma 6.8"`)。

原文 (Yoshida08 p.15):
> Lemma 6.8. For all σ ∈ G, we have i(σ[bar]) = 1|H| _τ∈H i(στ).

★`pdftotext` は分数の横線と総和記号 `Σ` を出力しない(構造化 HTML が
`data-txt="1|H|"` / `data-txt=""` で明示している)。読み下すと

  `i(σ̄) = (1/|H|) Σ_{τ∈H} i(στ)`   (`σ̄ = σH ∈ G/H`, `i(σ̄)` は `K''` で測る)。

設定は §6.2 (`#setup-6-2`):
> 6.2. The Hasse-Arf theorem. Let G = Gal(K′/K) with K′/K totally ramified as before,
> and let G ▷ H with G/H = Gal(K′′/K). For σ ∈ G, let σ[bar] = σH ∈ G/H be its image.

## ★原典の Proof(`.txt` 1033–1053 行)の 2 段が、そのまま §1 の 2 本になっている

> Let the minimal polynomial of π′ over O_{K′′} be f = ∏_{τ∈H}(X − τ(π′)) ∈ O_{K′′}[X].
> Applying σ, we get f^σ = ∏_{τ∈H}(X − στ(π′)) … Hence f^σ(π′) = ∏_{τ∈H}(π′ − στ(π′)) = ±b.
> First we prove a | b. As O_{K′′} = O[π′′], we have a | σ(x) − x for any x ∈ O_{K′′},
> hence a | f^σ − f, therefore a | f^σ(π′) − f(π′) = ±b. Now we prove b | a.
> Write π′′ = g(π′) for g ∈ O[X]. The polynomial g(X) − π′′ ∈ O_{K′′}[X] has π′ as a root,
> hence divisible by f in O_{K′′}[X]. Applying σ, we have f^σ | g(X) − σ(π′′) in O_{K′′}[X],
> hence g(π′) − σ(π′′) = −a is divisible by f^σ(π′) = ±b.

`a := σ(π′′) − π′′`、`b := ∏_{τ∈H}(στ(π′) − π′)` と置いて `a | b` と `b | a` を示す 2 段である。
★**この 2 段には分岐も付値も Galois 理論も出てこない**。出てくるのは
「多項式・環準同型・単項生成」だけなので、§1 で**一般の可換環**に切り出した。

## ★★除算を書かない —— 掛け算形で述べる

`ramIndex` は `ℕ∞` に値を取る。★`ℕ∞` に除算は無く、切り詰め引き算は「空虚に真」を作る
(`tools/lean-idioms.md` **#102**)。したがって本ファイルは原文の `1/|H|` を移項した

```
(Nat.card H : ℕ∞) * ramIndex ϖ σ = ∑ τ : H, ramIndex π' (σ * τ)
```

の形で述べる(`card_mul_ramIndex_eq_sum_ramIndex`)。★これは Y6 が Lemma 6.5 の
`−1` のずれに対して採った「`1 + v(…)` の形で持つ」処方と同じである。

## ★★商群 `G ⧸ H` を作っていない

原典は `σ̄ ∈ G/H` の `i(σ̄)` を `K''` 上で測る。本ファイルは
**`G` を `C`(= `O_{K′′}` に対応する DVR)に直接作用させ**、`ramIndex ϖ (σ : G)` と書く。
`H` が `C` に自明に作用する(`hHtriv`)ことを仮定に置いてあるので、この作用は自動的に
`G/H` を経由する。★**`MulSemiringAction (G ⧸ H) C` の構成は一切要らなかった。**
同じ理由で **`H` の正規性 `G ▷ H` も明示的には要らない** ——
正規性が要る内容(σ が `K′′` を保つこと)は `[MulSemiringAction G C]` と
`hcomp`(`algebraMap C B` が `G`-同変)に packaging されている。

## 本ファイルの構成

**§1 抽象核 —— 一般の可換環(分岐・付値・Galois の語彙が 1 つも出てこない)**

* `dvd_map_sub_self_of_mem_adjoin` : `φ` が `A` の像を固定するとき、`y ∈ A[t]` なら
  `φ t − t ∣ φ y − y`。★原文の「As O_{K′′} = O[π′′], we have a | σ(x) − x」。
* `dvd_eval_map_sub_eval` : 係数ごとに `a ∣ φ(c) − c` なら `a ∣ f^φ(z) − f(z)`。
  ★原文の「hence a | f^σ − f, therefore a | f^σ(π′) − f(π′)」。
* `map_prod_X_sub_C` / `eval_prod_X_sub_C` / `prod_sub_comm_eq` /
  `associated_prod_sub_comm` : `f^φ = ∏(X − φ(x i))` と `f^φ(z) = ±∏(φ(x i) − z)`(符号の `±`)。
* `dvd_prod_sub_map_of_coeff` : **`a ∣ b`**(第 1 段)。
* `prod_X_sub_C_dvd_of_forall_isRoot` : 相異なる根を全て持つ多項式は `∏(X − x i)` で割れる
  (整域)。★原文の「has π′ as a root, hence divisible by f」を**仮定に置かず導いた**。
* `prod_sub_map_dvd` : **`b ∣ a`**(第 2 段)。
* ★`associated_map_sub_self_prod` : **抽象核の本体**
  `Associated (φ w − w) (∏ i ∈ s, (φ (x i) − z))`。

**§2 具体層 —— 群作用と `Associated`**

* `associated_smul_sub_prod` : `Associated (σ • π'' − π'') (∏ τ : H, ((σ * τ) • π' − π'))`。

**§3 付値層 —— Lemma 6.8 そのもの**

* `addVal_smul_sub_eq_sum_ramIndex` : `v_{K′}(σ π'' − π'') = Σ_{τ∈H} i(στ)`。
* ★`card_mul_ramIndex_eq_sum_ramIndex` : **主結論**(掛け算形の Lemma 6.8)。

## ★★逸脱の記録

1. ★**`hram : ∀ c : C, addVal B (algebraMap C B c) = |H| * addVal C c` を仮定に置いた。**
   これは `v_{K′}|_{K′′} = |H| · v_{K′′}`、すなわち `K′/K′′` が完全分岐で `e(K′/K′′) = |H|`
   であることの言い換えである。★**この木にも mathlib にも在庫が無い**ことを
   `.cache/decl-index.txt` と `.cache/mathlib-index.txt` の両方で確認した
   (mathlib の `Ideal.ramificationIdx` は `addVal` と繋がっていない)。
   ★**導いていない。仮定である。** 下流(Proposition 6.9 Herbrand)が供給するか、
   別ノードとして立てる必要がある。
2. ★**`hfix : ∀ y : B, (∀ ρ ∈ H, ρ • y = y) → y ∈ Algebra.adjoin A {π''}` の 1 本に
   まとめた。** 原典はこれを「`B^H = O_{K′′}`」(Galois 理論)と
   「`O_{K′′} = O[π′′]`」(Lemma 5.11)の 2 つに分けて持っている。証明が使うのは
   合成の形だけなので 1 本にした。★弱めてはいない(2 つを仮定すれば従う)。
3. `hadj : Algebra.adjoin A {π'} = ⊤` は原典の `O_{K′} = O[π′]`(Lemma 5.11)。
   ★原典は「take uniformizers π′ and π′′ … so that O_{K′} = O[π′] and O_{K′′} = O[π′′]
   by Lemma 5.11」と書いており、仮定として明示されている。
   ★★本ファイルは `π'` が素元であることを**使っていない**(`hadj` しか使わない)ので、
   原典より仮定が弱い。
4. `C` を `B` の部分環として構成せず、`[Algebra C B]` と `hcomp` で受けた。
   `FixedPoints.subring` で `B^H` を作ると「その固定環が DVR である」ことの証明が
   別に要り、それは本補題の内容ではない。★`hϖ : algebraMap C B ϖ = π''` で
   両者を結び付けてある。

## 退化の自己検査

* ★`H` の有限性(`[Fintype ↥H]`)を落とすと `∏` も `∑` も書けない。
* ★`σ = id` のとき両辺 `⊤`(原文が "we understand the equality as ∞ = ∞" と断っている)。
  ★`ℕ∞` の掛け算は `0 * ⊤ = 0` なので、`Nat.card H = 0`(`H` が無限)だと
  **左辺が `0` に潰れて空虚に真**になる。`[Fintype ↥H]` から `0 < Nat.card H` が出るので
  この穴は塞がっている(`pos_natCard_subgroup`、`card_mul_ramIndex_one_eq_top`)。
* `hfix`(`C = A[π'']` 単項生成)を落とすと `a ∣ σx − x` が出ない(§1 第 1 段が崩れる)。
* `hram` を落とすと `|H|` の係数が根拠を失う(§3 が崩れる)。§1・§2 は影響を受けない。
-/

namespace ABC3.Found.PGC

open IsLocalRing ABC3.Skeleton.PGC IsDiscreteValuationRing
open Polynomial

def card_mul_ramIndex_eq_sum_ramIndex.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 15, item := "Lemma 6.8", sectionId := "lemma-6-8" }

/-! ## §1 抽象核 —— 一般の可換環

★分岐・付値・Galois の語彙は 1 つも現れない。原文の 2 段(`a | b` と `b | a`)が
そのまま 2 本の補題になる。 -/

/-- **原文「As O_{K′′} = O[π′′], we have a | σ(x) − x for any x ∈ O_{K′′}」** ——
環準同型 `φ` が `A` の像を各点固定するなら、`y ∈ A[t]` について `φ t − t ∣ φ y − y`。

★`y = P(t)`(`P` は `A` 係数)と書けば `φ y = P(φ t)` で、`P(u) − P(v)` は `u − v` で割れる
(`Polynomial.sub_dvd_eval_sub`)。 -/
theorem dvd_map_sub_self_of_mem_adjoin {A R : Type*} [CommRing A] [CommRing R] [Algebra A R]
    (φ : R →+* R) (hφ : ∀ a : A, φ (algebraMap A R a) = algebraMap A R a) {t y : R}
    (hy : y ∈ Algebra.adjoin A ({t} : Set R)) : φ t - t ∣ φ y - y := by
  rw [Algebra.adjoin_singleton_eq_range_aeval] at hy
  obtain ⟨q, hq⟩ := hy
  have hq' : (aeval t) q = y := hq
  set P : R[X] := q.map (algebraMap A R) with hP
  have hy1 : y = P.eval t := by rw [← hq', aeval_def, eval₂_eq_eval_map]
  have hPφ : P.map φ = P := by
    rw [hP, Polynomial.map_map]
    congr 1
    ext a
    simp [hφ]
  have h2 : φ (P.eval t) = P.eval (φ t) := by
    rw [← Polynomial.eval₂_at_apply (p := P) φ t, Polynomial.eval₂_eq_eval_map, hPφ]
  rw [hy1, h2]
  exact Polynomial.sub_dvd_eval_sub (φ t) t P

/-- **原文「hence a | f^σ − f, therefore a | f^σ(π′) − f(π′)」** ——
`f` の係数がすべて `a` を法として `φ` で不変なら、任意の点 `z` で
`a ∣ f^φ(z) − f(z)`。

★定数多項式による割り切りは係数ごとの割り切りと同値(`Polynomial.C_dvd_iff_dvd_coeff`)。 -/
theorem dvd_eval_map_sub_eval {R : Type*} [CommRing R] {a z : R} (φ : R →+* R) {f : R[X]}
    (hc : ∀ n, a ∣ φ (f.coeff n) - f.coeff n) :
    a ∣ (f.map φ).eval z - f.eval z := by
  have h : (C a) ∣ (f.map φ - f) := by
    rw [Polynomial.C_dvd_iff_dvd_coeff]
    intro i
    simpa [Polynomial.coeff_map] using hc i
  obtain ⟨q, hq⟩ := h
  have h3 := congrArg (Polynomial.eval z) hq
  simp only [Polynomial.eval_sub, Polynomial.eval_mul, Polynomial.eval_C] at h3
  exact ⟨q.eval z, h3⟩

/-- **原文「Applying σ, we get f^σ = ∏_{τ∈H}(X − στ(π′))」**。 -/
theorem map_prod_X_sub_C {R ι : Type*} [CommRing R] (s : Finset ι) (x : ι → R) (φ : R →+* R) :
    (∏ i ∈ s, (X - C (x i))).map φ = ∏ i ∈ s, (X - C (φ (x i))) := by
  rw [Polynomial.map_prod]
  simp

/-- `∏(X − x i)` の `z` での値は `∏(z − x i)`。 -/
theorem eval_prod_X_sub_C {R ι : Type*} [CommRing R] (s : Finset ι) (x : ι → R) (z : R) :
    (∏ i ∈ s, (X - C (x i))).eval z = ∏ i ∈ s, (z - x i) := by
  rw [Polynomial.eval_prod]
  simp

/-- 原文の `±` の中身 —— `∏(y i − z) = (−1)^{|s|} ∏(z − y i)`。 -/
theorem prod_sub_comm_eq {R ι : Type*} [CommRing R] (s : Finset ι) (y : ι → R) (z : R) :
    ∏ i ∈ s, (y i - z) = (-1) ^ s.card * ∏ i ∈ s, (z - y i) := by
  have h : ∀ i ∈ s, y i - z = (-1) * (z - y i) := fun i _ => by ring
  rw [Finset.prod_congr rfl h, Finset.prod_mul_distrib, Finset.prod_const]

/-- 原文の `f^σ(π′) = ±b` —— 符号は単位なので `Associated` で吸収する。 -/
theorem associated_prod_sub_comm {R ι : Type*} [CommRing R] (s : Finset ι) (y : ι → R) (z : R) :
    Associated (∏ i ∈ s, (z - y i)) (∏ i ∈ s, (y i - z)) :=
  ⟨(-1 : Rˣ) ^ s.card, by rw [prod_sub_comm_eq]; push_cast; ring⟩

/-- **第 1 段 `a | b`** —— `f := ∏(X − x i)` の係数がすべて `a` を法として `φ` で不変で、
`z` が `f` の根(`z = x i` なる `i ∈ s` がある)なら `a ∣ ∏(z − φ(x i))`。 -/
theorem dvd_prod_sub_map_of_coeff {R ι : Type*} [CommRing R] {s : Finset ι} {x : ι → R} {z a : R}
    (φ : R →+* R)
    (hc : ∀ n, a ∣ φ ((∏ i ∈ s, (X - C (x i))).coeff n) - (∏ i ∈ s, (X - C (x i))).coeff n)
    (hroot : ∃ i ∈ s, x i = z) :
    a ∣ ∏ i ∈ s, (z - φ (x i)) := by
  set f : R[X] := ∏ i ∈ s, (X - C (x i)) with hf
  have h1 : a ∣ (f.map φ).eval z - f.eval z := dvd_eval_map_sub_eval φ hc
  have h2 : f.eval z = 0 := by
    rw [hf, eval_prod_X_sub_C]
    obtain ⟨i, hi, hxi⟩ := hroot
    exact Finset.prod_eq_zero hi (by rw [hxi]; ring)
  have h3 : (f.map φ).eval z = ∏ i ∈ s, (z - φ (x i)) := by
    rw [hf, map_prod_X_sub_C, eval_prod_X_sub_C]
  rwa [h2, h3, sub_zero] at h1

/-- **第 2 段 `b | a`** —— `g^φ = g` かつ `w = g(z)` で `∏(X − x i) ∣ g − C w` なら
`∏(z − φ(x i)) ∣ φ w − w`。

★原文の「Applying σ, we have f^σ | g(X) − σ(π′′), hence g(π′) − σ(π′′) = −a is divisible by
f^σ(π′)」。`−a` の符号はここで吸収する。 -/
theorem prod_sub_map_dvd {R ι : Type*} [CommRing R] {s : Finset ι} {x : ι → R} {z w : R}
    (φ : R →+* R) {g : R[X]} (hg : g.map φ = g) (hw : w = g.eval z)
    (hdvd : (∏ i ∈ s, (X - C (x i))) ∣ g - C w) :
    (∏ i ∈ s, (z - φ (x i))) ∣ φ w - w := by
  have h1 : (∏ i ∈ s, (X - C (φ (x i)))) ∣ (g - C w).map φ := by
    rw [← map_prod_X_sub_C]
    exact Polynomial.map_dvd φ hdvd
  rw [Polynomial.map_sub, hg, Polynomial.map_C] at h1
  have h2 := Polynomial.eval_dvd (x := z) h1
  rw [eval_prod_X_sub_C, Polynomial.eval_sub, Polynomial.eval_C, ← hw] at h2
  obtain ⟨c, hc⟩ := h2
  exact ⟨-c, by rw [mul_neg, ← hc]; ring⟩

/-- **原文「has π′ as a root, hence divisible by f in O_{K′′}[X]」を仮定に置かずに導く。**

整域で `x` が `s` 上単射で、`q` がすべての `x i`(`i ∈ s`)を根に持つなら
`∏_{i∈s}(X − x i) ∣ q`。★`q = 0` の場合は `dvd_zero` で処理する
(`Multiset.prod_X_sub_C_dvd_iff_le_roots` は `q ≠ 0` を要求する)。 -/
theorem prod_X_sub_C_dvd_of_forall_isRoot {R ι : Type*} [CommRing R] [IsDomain R] {s : Finset ι}
    {x : ι → R} (hinj : ∀ i ∈ s, ∀ j ∈ s, x i = x j → i = j) {q : R[X]}
    (hq : ∀ i ∈ s, q.eval (x i) = 0) :
    (∏ i ∈ s, (X - C (x i))) ∣ q := by
  rcases eq_or_ne q 0 with rfl | hq0
  · exact dvd_zero _
  have hm : (Multiset.map (fun a => X - C a) (s.val.map x)).prod = ∏ i ∈ s, (X - C (x i)) := by
    rw [Multiset.map_map]
    rfl
  rw [← hm, Multiset.prod_X_sub_C_dvd_iff_le_roots hq0,
    Multiset.le_iff_subset (Multiset.Nodup.map_on hinj s.nodup)]
  intro a ha
  obtain ⟨i, hi, rfl⟩ := Multiset.mem_map.1 ha
  exact (Polynomial.mem_roots hq0).2 (hq i hi)

/-- ★★★**抽象核の本体(Yoshida 2008 Lemma 6.8 の中身)** ——

一般の可換整域 `R`、環準同型 `φ : R →+* R`、有限添字 `s`、`x : ι → R`、`z w : R` に対し

* `f := ∏_{i∈s}(X − x i)` の係数が `φ w − w` を法として `φ` で不変(`hc`)
* `z` は `f` の根(`hroot`)
* `g^φ = g`、`w = g(z)`、`f ∣ g − C w`(`hg`, `hw`, `hdvd`)

なら `Associated (φ w − w) (∏_{i∈s}(φ (x i) − z))`。

★**分岐・付値・Galois の語彙は 1 つも現れない。** 原典の設定(`a := σπ'' − π''`,
`b := ∏_{τ∈H}(στπ' − π')`)は `φ := σ`, `w := π''`, `x τ := τπ'`, `z := π'` を
代入するだけで得られる(§2)。 -/
theorem associated_map_sub_self_prod {R ι : Type*} [CommRing R] [IsDomain R] {s : Finset ι}
    {x : ι → R} {z w : R} (φ : R →+* R)
    (hc : ∀ n, (φ w - w) ∣ φ ((∏ i ∈ s, (X - C (x i))).coeff n)
      - (∏ i ∈ s, (X - C (x i))).coeff n)
    (hroot : ∃ i ∈ s, x i = z)
    {g : R[X]} (hg : g.map φ = g) (hw : w = g.eval z)
    (hdvd : (∏ i ∈ s, (X - C (x i))) ∣ g - C w) :
    Associated (φ w - w) (∏ i ∈ s, (φ (x i) - z)) :=
  (associated_of_dvd_dvd (dvd_prod_sub_map_of_coeff φ hc hroot)
    (prod_sub_map_dvd φ hg hw hdvd)).trans (associated_prod_sub_comm s (fun i => φ (x i)) z)

/-! ## §2 具体層 —— 群作用への代入

★§1 の核に `φ := σ`, `x τ := τ • π'`, `z := π'`, `w := π''` を代入する。 -/

/-- `G` が `B` に `A`-線型に作用するなら `σ` は `A` の像を各点固定する。 -/
theorem smul_algebraMap_self {A B : Type*} [CommRing A] [CommRing B] [Algebra A B]
    {G : Type*} [Group G] [MulSemiringAction G B] [SMulCommClass G A B] (σ : G) (a : A) :
    σ • algebraMap A B a = algebraMap A B a := by
  rw [Algebra.algebraMap_eq_smul_one, smul_comm, smul_one]

/-- `A` 係数の多項式は `σ` で不変。 -/
theorem map_algebraMap_toRingHom {A B : Type*} [CommRing A] [CommRing B] [Algebra A B]
    {G : Type*} [Group G] [MulSemiringAction G B] [SMulCommClass G A B] (σ : G) (g : A[X]) :
    (g.map (algebraMap A B)).map (MulSemiringAction.toRingHom G B σ) = g.map (algebraMap A B) := by
  rw [Polynomial.map_map]
  congr 1
  ext a
  exact smul_algebraMap_self σ a

/-- `A` 係数の多項式の値は `σ` と可換 —— `σ • g(y) = g(σ • y)`。 -/
theorem smul_eval_map {A B : Type*} [CommRing A] [CommRing B] [Algebra A B]
    {G : Type*} [Group G] [MulSemiringAction G B] [SMulCommClass G A B] (σ : G) (g : A[X])
    (y : B) :
    σ • ((g.map (algebraMap A B)).eval y) = (g.map (algebraMap A B)).eval (σ • y) := by
  have h := Polynomial.eval₂_at_apply (p := g.map (algebraMap A B))
    (MulSemiringAction.toRingHom G B σ) y
  rw [Polynomial.eval₂_eq_eval_map, map_algebraMap_toRingHom] at h
  exact h.symm

/-- **原文「f = ∏_{τ∈H}(X − τ(π′)) ∈ O_{K′′}[X]」の "∈ O_{K′′}[X]" の部分** ——
`ρ ∈ H` は `H` を平行移動で置換するので `f` を保つ。 -/
theorem map_conjProd_self {B : Type*} [CommRing B] {G : Type*} [Group G] [MulSemiringAction G B]
    {H : Subgroup G} [Fintype H] {π' : B} (ρ : G) (hρ : ρ ∈ H) :
    (∏ τ : H, (X - C ((τ : G) • π'))).map (MulSemiringAction.toRingHom G B ρ)
      = ∏ τ : H, (X - C ((τ : G) • π')) := by
  rw [map_prod_X_sub_C]
  have hstep : ∀ τ : H, (MulSemiringAction.toRingHom G B ρ) ((τ : G) • π')
      = (((Equiv.mulLeft (⟨ρ, hρ⟩ : H)) τ : H) : G) • π' := by
    intro τ
    show ρ • ((τ : G) • π') = _
    rw [← mul_smul]
    rfl
  calc ∏ τ : H, (X - C ((MulSemiringAction.toRingHom G B ρ) ((τ : G) • π')))
      = ∏ τ : H, (X - C ((((Equiv.mulLeft (⟨ρ, hρ⟩ : H)) τ : H) : G) • π')) :=
        Finset.prod_congr rfl fun τ _ => by rw [hstep τ]
    _ = ∏ τ : H, (X - C ((τ : G) • π')) :=
        Equiv.prod_comp (Equiv.mulLeft (⟨ρ, hρ⟩ : H)) (fun τ : H => X - C ((τ : G) • π'))

/-- `f` の各係数は `H` で不変(したがって `B^H` に属する)。 -/
theorem smul_coeff_conjProd {B : Type*} [CommRing B] {G : Type*} [Group G] [MulSemiringAction G B]
    {H : Subgroup G} [Fintype H] {π' : B} (ρ : G) (hρ : ρ ∈ H) (n : ℕ) :
    ρ • ((∏ τ : H, (X - C ((τ : G) • π'))).coeff n)
      = (∏ τ : H, (X - C ((τ : G) • π'))).coeff n := by
  have h2 : (MulSemiringAction.toRingHom G B ρ) ((∏ τ : H, (X - C ((τ : G) • π'))).coeff n)
      = (∏ τ : H, (X - C ((τ : G) • π'))).coeff n := by
    rw [← Polynomial.coeff_map, map_conjProd_self ρ hρ]
  exact h2

/-- `τ ↦ τ • π'` は `H` 上単射 —— `π'` が `B` を `A` 上生成することから出る。
★§1 の `prod_X_sub_C_dvd_of_forall_isRoot`(根の相異性)に渡すために要る。 -/
theorem coe_smul_injOn {A B : Type*} [CommRing A] [CommRing B] [Algebra A B]
    {G : Type*} [Group G] [MulSemiringAction G B] [SMulCommClass G A B] [FaithfulSMul G B]
    {H : Subgroup G} {π' : B} (hadj : Algebra.adjoin A ({π'} : Set B) = ⊤) {τ₁ τ₂ : H}
    (h : (τ₁ : G) • π' = (τ₂ : G) • π') : τ₁ = τ₂ := by
  have h1 : ((τ₂ : G)⁻¹ * (τ₁ : G)) • π' = π' := by
    rw [mul_smul, h, inv_smul_smul]
  have h2 : ((τ₂ : G)⁻¹ * (τ₁ : G)) = 1 :=
    eq_one_of_smul_eq_of_adjoin_eq_top (A := A) hadj h1
  exact Subtype.ext (inv_mul_eq_one.mp h2).symm

/-- ★★**具体層の本体** —— `a := σ • π'' − π''` と `b := ∏_{τ∈H}((στ) • π' − π')` は
同伴である(原文の `v_{K′}(a) = v_{K′}(b)`、割り切りの両向き)。

* `hadj` : `O_{K′} = O[π′]`(原典 Lemma 5.11)
* `hfix` : `H` で固定される元はすべて `A[π'']` に属する
  (原典の `B^H = O_{K′′}` と `O_{K′′} = O[π′′]` を合わせたもの)
* `hinv` : `π''` は `H` で固定される

★`π'`・`π''` が素元であることは使っていない。 -/
theorem associated_smul_sub_prod {A B : Type*} [CommRing A] [CommRing B] [Algebra A B]
    [IsDomain B] {G : Type*} [Group G] [MulSemiringAction G B] [SMulCommClass G A B]
    [FaithfulSMul G B] {H : Subgroup G} [Fintype H] {π' π'' : B}
    (hadj : Algebra.adjoin A ({π'} : Set B) = ⊤)
    (hfix : ∀ y : B, (∀ ρ ∈ H, ρ • y = y) → y ∈ Algebra.adjoin A ({π''} : Set B))
    (hinv : ∀ ρ ∈ H, ρ • π'' = π'') (σ : G) :
    Associated (σ • π'' - π'') (∏ τ : H, ((σ * (τ : G)) • π' - π')) := by
  classical
  obtain ⟨g, hg⟩ : ∃ g : A[X], (aeval π' g : B) = π'' := by
    have hmem : π'' ∈ Algebra.adjoin A ({π'} : Set B) := by rw [hadj]; trivial
    rw [Algebra.adjoin_singleton_eq_range_aeval] at hmem
    obtain ⟨g, hg⟩ := hmem
    exact ⟨g, hg⟩
  have hPeval : (g.map (algebraMap A B)).eval π' = π'' := by
    rw [← eval₂_eq_eval_map, ← aeval_def, hg]
  have hgmap : (g.map (algebraMap A B)).map (MulSemiringAction.toRingHom G B σ)
      = g.map (algebraMap A B) := map_algebraMap_toRingHom σ g
  have hdvd : (∏ τ : H, (X - C ((τ : G) • π'))) ∣ (g.map (algebraMap A B)) - C π'' := by
    refine prod_X_sub_C_dvd_of_forall_isRoot (fun i _ j _ h => coe_smul_injOn (A := A) hadj h) ?_
    intro τ _
    have hev : (g.map (algebraMap A B)).eval ((τ : G) • π') = π'' := by
      rw [← smul_eval_map (τ : G) g π', hPeval]
      exact hinv (τ : G) τ.2
    simp [hev]
  have hc : ∀ n, ((MulSemiringAction.toRingHom G B σ) π'' - π'')
      ∣ (MulSemiringAction.toRingHom G B σ) ((∏ τ : H, (X - C ((τ : G) • π'))).coeff n)
        - (∏ τ : H, (X - C ((τ : G) • π'))).coeff n := by
    intro n
    refine dvd_map_sub_self_of_mem_adjoin (A := A) _ (fun a => smul_algebraMap_self σ a) ?_
    exact hfix _ (fun ρ hρ => smul_coeff_conjProd ρ hρ n)
  have hroot : ∃ τ ∈ (Finset.univ : Finset H), (τ : G) • π' = π' :=
    ⟨1, Finset.mem_univ _, by simp⟩
  have hmain := associated_map_sub_self_prod (MulSemiringAction.toRingHom G B σ) hc hroot
    hgmap hPeval.symm hdvd
  have hrw : (∏ τ : H, ((MulSemiringAction.toRingHom G B σ) ((τ : G) • π') - π'))
      = ∏ τ : H, ((σ * (τ : G)) • π' - π') :=
    Finset.prod_congr rfl fun τ _ => by
      show σ • ((τ : G) • π') - π' = _
      rw [mul_smul]
  rw [hrw] at hmain
  exact hmain

/-! ## §3 付値層 —— Lemma 6.8

★`ℕ∞` の除算・引き算は使わない(掛け算形)。 -/

/-- 有限積の付値は付値の和(`addVal` の `map_mul` を `Finset` に持ち上げただけ)。 -/
theorem addVal_prod {B : Type*} [CommRing B] [IsDomain B] [IsDiscreteValuationRing B]
    {ι : Type*} (s : Finset ι) (f : ι → B) :
    addVal B (∏ i ∈ s, f i) = ∑ i ∈ s, addVal B (f i) := by
  classical
  induction s using Finset.induction_on with
  | empty => simp
  | @insert a s ha ih => rw [Finset.prod_insert ha, Finset.sum_insert ha, addVal_mul, ih]

/-- **原文の `v_{K′}(a) = v_{K′}(b)`** —— `v_{K′}(σ π'' − π'') = Σ_{τ∈H} i(στ)`。

★これが Lemma 6.8 の内容そのもの(`K′` 側だけで書いた形)。`|H|` の因子は
§3 の最後で `v_{K′}|_{K′′} = |H| v_{K′′}` を使って現れる。 -/
theorem addVal_smul_sub_eq_sum_ramIndex {A B : Type*} [CommRing A] [CommRing B] [Algebra A B]
    [IsDomain B] [IsDiscreteValuationRing B] {G : Type*} [Group G] [MulSemiringAction G B]
    [SMulCommClass G A B] [FaithfulSMul G B] {H : Subgroup G} [Fintype H] {π' π'' : B}
    (hadj : Algebra.adjoin A ({π'} : Set B) = ⊤)
    (hfix : ∀ y : B, (∀ ρ ∈ H, ρ • y = y) → y ∈ Algebra.adjoin A ({π''} : Set B))
    (hinv : ∀ ρ ∈ H, ρ • π'' = π'') (σ : G) :
    addVal B (σ • π'' - π'') = ∑ τ : H, ramIndex π' (σ * (τ : G)) := by
  rw [(addVal_eq_iff_associated _ _).2 (associated_smul_sub_prod (A := A) hadj hfix hinv σ),
    addVal_prod]
  rfl

/-- ★退化の自己検査 —— `H` が有限なら `0 < |H|`。これが無いと `ℕ∞` の
`0 * ⊤ = 0` で主結論の左辺が潰れ、空虚に真になる。 -/
theorem pos_natCard_subgroup {G : Type*} [Group G] (H : Subgroup G) [Fintype H] :
    0 < Nat.card H := Nat.card_pos

/-- ★退化の自己検査 —— `σ = id` では左辺が `⊤`(原文の "we understand the equality as
∞ = ∞")。`0 < |H|` があるので `|H| * ⊤ = ⊤` になる。 -/
theorem card_mul_ramIndex_one_eq_top {C : Type*} [CommRing C] [IsDomain C]
    [IsDiscreteValuationRing C] {G : Type*} [Group G] [MulSemiringAction G C]
    (H : Subgroup G) [Fintype H] (ϖ : C) :
    (Nat.card H : ℕ∞) * ramIndex ϖ (1 : G) = ⊤ := by
  rw [ramIndex_one]
  exact ENat.mul_top (by exact_mod_cast (Nat.card_pos (α := H)).ne')

/-- ★★★★**Yoshida 2008 Lemma 6.8**(掛け算形) ——

`|H| · i(σ̄) = Σ_{τ∈H} i(στ)`。

原文 (Yoshida08 p.15):
> Lemma 6.8. For all σ ∈ G, we have i(σ[bar]) = 1|H| _τ∈H i(στ).

* `C` は `O_{K′′}` に対応する DVR。`G` が `C` に作用し、`H` は自明に作用する(`hHtriv`)ので
  この作用は自動的に `G/H` を経由する。★**商群 `G ⧸ H` を作る必要は無い。**
* `hcomp` : `algebraMap C B` は `G`-同変(`K′′ ⊂ K′` が `G`-安定であること)。
* `hram` : ★**`v_{K′}|_{K′′} = |H| · v_{K′′}`(`K′/K′′` が完全分岐で `e = |H|`)。
  これは在庫に無いので仮定に置いた。導いていない。**
* `hϖ` : `ϖ ∈ C` が `π'' ∈ B` に写る(`K′′` の素元)。
* `hadj` / `hfix` : 原典 Lemma 5.11 の `O_{K′} = O[π′]` / `O_{K′′} = O[π′′]`。

★原文の `1/|H|` は書かない。`ℕ∞` に除算は無く、切り詰め引き算は空虚な真を作る
(`lean-idioms.md` #102)。 -/
theorem card_mul_ramIndex_eq_sum_ramIndex {A B : Type*} [CommRing A] [CommRing B] [Algebra A B]
    [IsDomain B] [IsDiscreteValuationRing B] {G : Type*} [Group G] [MulSemiringAction G B]
    [SMulCommClass G A B] [FaithfulSMul G B] {H : Subgroup G} [Fintype H] {π' π'' : B}
    {C : Type*} [CommRing C] [IsDomain C] [IsDiscreteValuationRing C] [Algebra C B]
    [MulSemiringAction G C]
    (hcomp : ∀ (ρ : G) (c : C), algebraMap C B (ρ • c) = ρ • algebraMap C B c)
    (hHtriv : ∀ ρ ∈ H, ∀ c : C, ρ • c = c)
    (hram : ∀ c : C, addVal B (algebraMap C B c) = (Nat.card H : ℕ∞) * addVal C c)
    {ϖ : C} (hϖ : algebraMap C B ϖ = π'')
    (hadj : Algebra.adjoin A ({π'} : Set B) = ⊤)
    (hfix : ∀ y : B, (∀ ρ ∈ H, ρ • y = y) → y ∈ Algebra.adjoin A ({π''} : Set B))
    (σ : G) :
    (Nat.card H : ℕ∞) * ramIndex ϖ σ = ∑ τ : H, ramIndex π' (σ * (τ : G)) := by
  have hinv : ∀ ρ ∈ H, ρ • π'' = π'' := by
    intro ρ hρ
    rw [← hϖ, ← hcomp ρ ϖ, hHtriv ρ hρ ϖ]
  have hkey : algebraMap C B (σ • ϖ - ϖ) = σ • π'' - π'' := by
    rw [map_sub, hcomp σ ϖ, hϖ]
  calc (Nat.card H : ℕ∞) * ramIndex ϖ σ
      = addVal B (algebraMap C B (σ • ϖ - ϖ)) := (hram _).symm
    _ = addVal B (σ • π'' - π'') := by rw [hkey]
    _ = ∑ τ : H, ramIndex π' (σ * (τ : G)) :=
        addVal_smul_sub_eq_sum_ramIndex (A := A) hadj hfix hinv σ

end ABC3.Found.PGC
