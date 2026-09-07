import ABC3.Found.PGC.LubinTateEndoTwisted

/-!
# 半線型方程式の完備版と、ねじれ Lubin-Tate 舞台の `K̂^ur` への当てはめ

典拠: T. Yoshida, *Local Class Field Theory via Lubin-Tate Theory*
(Ann. Fac. Sci. Toulouse 17-2, 2008; arXiv math/0606108) Lemma 3.4(物理 p.5)。
構造化済み原文は
`ResearchPaper/1_Structured/Local Class Field Theory via Lubin-Tate Theory/`
の `section-3.html` の `#lemma-3-4`。

原文 (Yoshida08 p.5):
> Lemma 3.4. Let π be a uniformizer of L, and let f ∈ O[scr]_L[[X]] satisfy: (3.2.1) f(X) ≡ πX (mod deg 2), f(X) ≡ X^q (mod p[frak]). Let π′, f′ be another such pair. Assume that θ_1, . . . , θ_t ∈ Θ^L_π,π′. Then there is a unique F ∈ O[scr]_L[[X_1, . . . , X_t]] satisfying the following: F ≡ θ_1X_1 + · · · + θ_tX_t (mod deg 2), f′ ◦ F = F^ϕ ◦ f.

## 直前の波との関係(★本ファイルの存在理由)

`Found/PGC/LubinTateEndoTwisted.lean` が `𝒪_L` 係数の `[θ]_{f,f′}` と
ねじれ版の関数等式 `f′ ∘ [θ] = [θ]^ϕ ∘ f` を閉じたが、その逸脱記録の第 1 項に

> 半線型方程式 `c − u·ϕ(c) = b` の可解性を舞台の仮定 `hsolve` に括り出した。
> 十分条件は「有限性 + 中山」と「`ϕ = id`」の 2 つで、
> 完備(非有限)な `L = K̂^ur` 用の `hsolve` は本ファイルに無い

と残っていた。★本ファイルはその 1 点を埋め、あわせて舞台を木の実物
(`unramifiedCompletionInt` / `arithFrobenius`)へ当てはめる。

## 原典のどの式か(★`.txt` 236–253 行の直読で確定した)

Lemma 3.4 の証明は次数 `m` の帰納で、`m+1` 次の単項式ごとに未知係数 `α` について

```
π′β + π′α − π^{m+1} α^ϕ = 0,   すなわち   α − (π^{m+1}/π′)·α^ϕ = −β
```

を解き、答を無限級数

```
α = −β − Σ_{i≥1} (π^{m+1}/π′)^{1+ϕ+⋯+ϕ^{i−1}} β^{ϕ^i}
```

で与える。★`u := π^{m+1}/π′` は `m ≥ 1` のとき付値 `m ≥ 1` すなわち `u ∈ 𝔪` であり、
半線型作用素 `T c := u·ϕ(c)` は**縮小**である。原典の無限級数はちょうど
`Σ_{i≥0} T^i(−β)` であり(`T^i(b) = u^{1+ϕ+⋯+ϕ^{i−1}}·b^{ϕ^i}`)、
その収束が `L` の完備性から出る、というのが原典の「Appendix I の収束」である。

★本ファイルはこの級数を **`IsAdicComplete` だけで**扱う。ノルムも距離も
`Tendsto` も現れない —— `IsPrecomplete.prec` で極限を取り、`IsHausdorff.haus` で
方程式を確定させる(Λ6 の `DworkAdditive` と同じ流儀)。

## 何を足したか

### 純抽象核(分岐・付値・Galois・Lubin-Tate のどれも出てこない)

| 宣言 | 内容 |
|---|---|
| `map_mem_pow_of_map_mem` | `ϕ(I) ⊆ I` ならば `ϕ(I^n) ⊆ I^n` |
| `semilinearGeomSum` | `S 0 = 0`, `S (n+1) = b + u·ϕ(S n)`(= `Σ_{i<n} T^i b`) |
| `semilinearGeomSum_sub_mem` | Cauchy 評価 `S(m+k) − S m ∈ I^m` |
| `exists_sub_mul_map_eq_of_isAdicComplete` | ★★**核 1**: `IsAdicComplete I A` と `ϕ(I) ⊆ I` と `u ∈ I` から `c − u·ϕ(c) = b` を解く |
| `semilinear_solution_unique_of_isHausdorff` | ★★**核 2**: `IsHausdorff I A` だけで解は高々 1 つ |
| `maximalIdeal_map_mem_of_ringEquiv` | ★**核 3**: 環同型は局所環の極大イデアルを保つ |

★核 1・核 2 は `CommRing A` と `Ideal A` と `A →+* A` しか使わない。
★★核 1 の再帰 `S (n+1) = b + u·ϕ(S n)` は**定義そのもの**なので、
`Σ_{i<n} T^i b` という `Finset` の和を一度も書かずに済んでいる
(原典の `1+ϕ+⋯+ϕ^{i−1}` という指数は Lean には現れない)。

### 配管

| 宣言 | 内容 |
|---|---|
| `hsolve_of_isAdicComplete` | 核 1 を `I = 𝔪` に固定して `TwistedLT.hsolve` の形に整えたもの |
| `hsolve_of_isAdicComplete_ringEquiv` | `ϕ` が環同型のとき `ϕ(𝔪) ⊆ 𝔪` を自動で埋める形 |
| `hsolve_unramifiedCompletionInt` | ★★**`𝒪_{K̂^ur}` での `hsolve`**(`σ` は任意の `unramGal K`) |
| `hsolve_arithFrobenius_unramifiedCompletionInt` | 同・算術 Frobenius 版 |
| `TwistedLT.ofAdicComplete` | 舞台の構成子(`hsolve` を完備性から自動で埋める) |
| `TwistedLT.ofUnramifiedCompletion` | 舞台を `A = 𝒪_{K̂^ur}`, `ϕ = σ` に当てはめたもの |
| `nonempty_twistedLT_unramifiedCompletionInt` | ★★**舞台が空でない**(`f = π′X + X^q` で実際に 1 つ作る) |

## `ϕ(𝔪) ⊆ 𝔪` をどこから出したか(★退化の自己検査)

核 1 は `ϕ(I) ⊆ I` を**仮定**として受け取る。木の実物ではこれを

* `unramGalCompletionInt K σ` が**環同型**(`≃+*`)であること

だけから出している(`maximalIdeal_map_mem_of_ringEquiv`)。すなわち
`IsLocalHom` の instance `isLocalHom_equiv`(`Mathlib/Algebra/Group/Units/Equiv.lean`)
により `IsUnit (e x) → IsUnit x` が成り立ち、`𝔪 = nonunits` の像は `nonunits` に入る。
★**付値も不分岐性も使っていない**——Galois 作用が全単射であることだけが効く。
★`σπ = π`(`unramGalCompletionInt_uniformizerCompletionInt`)は**要らなかった**
(Λ6 の加法版 Dwork はそれを使うが、こちらは使わない)。

## ★★加法版 Dwork はそのまま当たるか —— 当たらない(実測)

`DworkAdditive.exists_unramGalCompletionInt_sub_self_eq` は

```
σb − b = c        すなわち   (ϕ − 1)b = c
```

を解く定理である。本ファイルが解くのは

```
c − u·ϕ(c) = b    すなわち   (1 − u·ϕ)c = b   (u ∈ 𝔪)
```

で、**別の方程式**である。違いは本質的で、向きが逆になっている:

* Dwork の `ϕ − 1` は**縮小ではない**。だから剰余体で `t^q − t = c̄` を解く段
  (`exists_pow_card_sub_self_eq_completion`、`𝓀_{K̂^ur}` が代数閉であること)が要る。
  ★不分岐性・剰余体の代数閉性という**分岐論の入力**が本質的に効く。
* 本ファイルの `1 − u·ϕ` は `u ∈ 𝔪` なので**縮小**であり、幾何級数がそのまま収束する。
  ★剰余体を一度も見ない。**純粋な可換環論**で閉じる。

★したがって Dwork は当てられなかった。★**当たらなかったのは「弱いから」ではなく
「Dwork の方が難しい問題を解いているから」**である。
一方で Dwork が作った道具は 3 つ効いた:

* `isAdicComplete_unramifiedCompletionInt`(`𝒪_{K̂^ur}` の `𝔪` 進完備性)——★これが無いと配管が立たない、
* `sModEq_pow_iff`(`SModEq` とイデアル所属の往復)、
* `residue_unramGalCompletionInt`(完備化側での `hϕres`)。

## 退化の自己検査

* ★★**`hsolve` の完備版は空虚でない**。`hsolve_unramifiedCompletionInt` は
  **仮定なしで** `A = 𝒪_{K̂^ur}`, `ϕ = unramGalCompletionInt K σ` に対し成り立つ
  (`σ` は任意の `unramGal K` でよい——ここでも剰余体の条件は要らない)。
  さらに `nonempty_twistedLT_unramifiedCompletionInt` で
  **`TwistedLT ↥(unramifiedCompletionInt K) p (absoluteInertiaDegree K)` の元を 1 つ作っている**
  (`f = f′ = π′X + X^{q}`、`q = Nat.card 𝓀[K]`)。★舞台全体が空でない。
* ★★**既存の 2 つの十分条件と整合する**。`semilinear_solution_subsingleton` が
  「Hausdorff なら解は高々 1 つ」を与えるので、`hsolve` の答えは十分条件の
  取り方によらない。とくに `eq_of_one_sub_mul_eq_of_isHausdorff` は、`ϕ = id` のとき
  完備版の解が `TwistedLT.ofUntwisted` の解 `(1−u)⁻¹b` と**一致する**ことを言う。
  有限性 + 中山版(`exists_semilinear_solution_of_moduleFinite`)についても、
  舞台が Hausdorff でありさえすれば同じ `c` になる(同じ核 2 の帰結)。
  ★重なる範囲で 2 つの答えが食い違わないことを、存在ではなく**一意性**で示した。
* ★`u ∈ 𝔪` を落とすと**偽**。`u = 1`, `ϕ = id` なら `c − c = 0 ≠ b`。
  縮小性がどこで効いているかは `semilinearGeomSum_sub_mem` の `pow_succ'` の 1 行である。
* ★`IsHausdorff` を落とすと一意性は**偽**(`I` 進位相が分離でないと `d = u·ϕ(d)` の
  非自明解がありうる)。`IsPrecomplete` を落とすと存在が**偽**。
  ★両方を使っており、`IsAdicComplete` のどちらの半分も遊んでいない。
* ★`ℕ∞` の切り詰め引き算・除算は 1 度も書いていない。原典の `π^{m+1}/π′` という
  **商**は、`u` を最初から `𝔪` の元として受け取ることで回避した(`#102`)。

## 逸脱の記録

1. ★原典は `u = π^{m+1}/π′` という**具体的な商**で書くが、本ファイルは
   `u ∈ 𝔪` という**性質だけ**を仮定する(一般化)。可除性
   (`π^{m+1} = π′u`)は `LubinTateEndoTwisted.exists_next_step_twisted` 側で
   既に処理済みなので、後続への影響は無い。
2. ★原典は `L` を「`K` の完備不分岐拡大」とし、収束を Appendix I の付値評価で
   議論する。本ファイルは付値を使わず `IsAdicComplete 𝔪 A` で置き換えた。
   `A = 𝒪_{K̂^ur}` では両者は一致する(`isAdicComplete_unramifiedCompletionInt`)。
   ★得られる主張は同じで、より一般の環にも当たる。
3. ★`TwistedLT.ofUnramifiedCompletion` は `f`・`f′` を**引数として受け取る**
   (原典のように「(3.2.1) を満たす `f` を 1 つ取る」とは書いていない)。
   実際に取れることは `nonempty_twistedLT_unramifiedCompletionInt` で別に示した。
4. Lemma 3.4 の `t` 変数版(`F ∈ 𝒪_L[[X_1,…,X_t]]`)は入っていない。
   本ファイルが埋めるのは証明中の**係数の方程式 1 本**であり、
   1 変数版の帰納は `LubinTateEndoTwisted` 側に既にある。

## 新しく必要になったノード

* Corollary 4.9 具体化の項目 2(完備化レベルの半線型性
  `σ([θ](α)) = [θ^{(j)}](σα)`)——点で評価する層。★本ファイルには入っていない。
* Proposition 3.5 の (i)(形式群 `F_f` の存在)・(iii)・加法性のねじれ版
  (2 変数の一意性補題のねじれ版が要る)。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC ABC3.Found.GaloisRep
open scoped NNReal Valued

def exists_sub_mul_map_eq_of_isAdicComplete.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 5, item := "Lemma 3.4", sectionId := "lemma-3-4" }

def semilinear_solution_unique_of_isHausdorff.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 5, item := "Lemma 3.4", sectionId := "lemma-3-4" }

def hsolve_of_isAdicComplete.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 5, item := "Lemma 3.4", sectionId := "lemma-3-4" }

def hsolve_unramifiedCompletionInt.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 5, item := "Lemma 3.4", sectionId := "lemma-3-4" }

def TwistedLT.ofAdicComplete.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 5, item := "Lemma 3.4", sectionId := "lemma-3-4" }

def nonempty_twistedLT_unramifiedCompletionInt.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 5, item := "Proposition 3.5", sectionId := "prop-3-5" }

/-! ## 0. 純抽象核

★この節には分岐・付値・Galois・Lubin-Tate・冪級数のいずれも出てこない。
一般の可換環 `A` とイデアル `I` と環準同型 `ϕ : A →+* A` しか使わない。 -/

/-- `ϕ(I) ⊆ I` ならば `ϕ(I^n) ⊆ I^n`。

★`Ideal.map_le_iff_le_comap` は項で書くと `?m` が決まらないので
`rw` してから `intro` する(`lean-idioms.md` #190)。 -/
theorem map_mem_pow_of_map_mem {A : Type*} [CommRing A] {I : Ideal A} (ϕ : A →+* A)
    (hϕ : ∀ x ∈ I, ϕ x ∈ I) (n : ℕ) {x : A} (hx : x ∈ I ^ n) : ϕ x ∈ I ^ n := by
  have hmap : Ideal.map ϕ I ≤ I := by
    rw [Ideal.map_le_iff_le_comap]
    intro y hy
    exact hϕ y hy
  have hpow : Ideal.map ϕ (I ^ n) ≤ I ^ n := by
    rw [Ideal.map_pow]
    exact Ideal.pow_right_mono hmap n
  exact hpow (Ideal.mem_map_of_mem ϕ hx)

/-- 半線型作用素 `T c := u·ϕ(c)` の幾何級数の部分和 `Σ_{i<n} T^i b`。

★`Finset` の和ではなく漸化式で定義する。すると原典の
`α = −β − Σ_{i≥1}(π^{m+1}/π′)^{1+ϕ+⋯+ϕ^{i−1}}β^{ϕ^i}` に現れる
**指数 `1+ϕ+⋯+ϕ^{i−1}` が Lean には一度も現れない**——
`S (n+1) = b + u·ϕ(S n)` が定義そのものだからである。 -/
def semilinearGeomSum {A : Type*} [CommRing A] (ϕ : A →+* A) (u b : A) : ℕ → A
  | 0 => 0
  | n + 1 => b + u * ϕ (semilinearGeomSum ϕ u b n)

@[simp] theorem semilinearGeomSum_zero {A : Type*} [CommRing A] (ϕ : A →+* A) (u b : A) :
    semilinearGeomSum ϕ u b 0 = 0 := rfl

@[simp] theorem semilinearGeomSum_succ {A : Type*} [CommRing A] (ϕ : A →+* A) (u b : A) (n : ℕ) :
    semilinearGeomSum ϕ u b (n + 1) = b + u * ϕ (semilinearGeomSum ϕ u b n) := rfl

/-- **Cauchy 評価**: 部分和は `I^m` を法として安定する。

★縮小性(`u ∈ I`)が効くのはここ 1 箇所だけである(`pow_succ'` の行)。 -/
theorem semilinearGeomSum_sub_mem {A : Type*} [CommRing A] {I : Ideal A} (ϕ : A →+* A)
    (hϕ : ∀ x ∈ I, ϕ x ∈ I) {u : A} (hu : u ∈ I) (b : A) (m k : ℕ) :
    semilinearGeomSum ϕ u b (m + k) - semilinearGeomSum ϕ u b m ∈ I ^ m := by
  induction m generalizing k with
  | zero => simp
  | succ m ih =>
      have hstep : semilinearGeomSum ϕ u b (m + 1 + k) - semilinearGeomSum ϕ u b (m + 1)
          = u * ϕ (semilinearGeomSum ϕ u b (m + k) - semilinearGeomSum ϕ u b m) := by
        have hidx : m + 1 + k = (m + k) + 1 := by omega
        rw [hidx, semilinearGeomSum_succ, semilinearGeomSum_succ, map_sub, mul_sub]
        ring
      rw [hstep, pow_succ' I m]
      exact Ideal.mul_mem_mul hu (map_mem_pow_of_map_mem ϕ hϕ m (ih k))

/-- ★★★★★★★★**純抽象核 1(完備性から半線型方程式を解く)**。

`I` 進完備な可換環 `A` の上で、`ϕ(I) ⊆ I` なる `ϕ : A →+* A` と `u ∈ I` に対し
半線型方程式 `c − u·ϕ(c) = b` は解ける。

★これが Yoshida08 Lemma 3.4 の証明で「Appendix I の収束」と呼ばれている段の中身
であり、ねじれ版の再帰(`LubinTateEndoTwisted.exists_next_step_twisted`)が要求する
唯一の新しい入力である。★分岐も付値も Lubin-Tate も出てこない。

証明は 3 行:`semilinearGeomSum` が Cauchy 列であること(`IsPrecomplete`)、
極限 `c` を取ること、`c − u·ϕ(c) − b ∈ I^n` を全ての `n` で示して `IsHausdorff`。 -/
theorem exists_sub_mul_map_eq_of_isAdicComplete {A : Type*} [CommRing A] (I : Ideal A)
    [IsAdicComplete I A] (ϕ : A →+* A) (hϕ : ∀ x ∈ I, ϕ x ∈ I) {u : A} (hu : u ∈ I) (b : A) :
    ∃ c : A, c - u * ϕ c = b := by
  obtain ⟨c, hc⟩ := IsPrecomplete.prec (I := I) (M := A) ‹IsAdicComplete I A›.toIsPrecomplete
    (f := semilinearGeomSum ϕ u b) (by
      intro m n hmn
      rw [sModEq_pow_iff]
      obtain ⟨k, rfl⟩ := Nat.exists_eq_add_of_le hmn
      simpa using Submodule.neg_mem _ (semilinearGeomSum_sub_mem ϕ hϕ hu b m k))
  refine ⟨c, ?_⟩
  have key : ∀ n : ℕ, c - u * ϕ c - b ∈ I ^ n := by
    intro n
    have h1 : c - semilinearGeomSum ϕ u b (n + 1) ∈ I ^ n := by
      have hmem := (sModEq_pow_iff I (n + 1) (semilinearGeomSum ϕ u b (n + 1)) c).mp (hc (n + 1))
      exact Ideal.pow_le_pow_right (Nat.le_succ n) (by simpa using Submodule.neg_mem _ hmem)
    have h2 : u * ϕ (semilinearGeomSum ϕ u b n - c) ∈ I ^ n := by
      have hmem := (sModEq_pow_iff I n (semilinearGeomSum ϕ u b n) c).mp (hc n)
      refine Ideal.pow_le_pow_right (Nat.le_succ n) ?_
      rw [pow_succ' I n]
      exact Ideal.mul_mem_mul hu (map_mem_pow_of_map_mem ϕ hϕ n hmem)
    have heq : c - u * ϕ c - b
        = (c - semilinearGeomSum ϕ u b (n + 1)) + u * ϕ (semilinearGeomSum ϕ u b n - c) := by
      rw [semilinearGeomSum_succ, map_sub, mul_sub]
      ring
    rw [heq]
    exact Ideal.add_mem _ h1 h2
  exact (IsHausdorff.eq_iff_smodEq (I := I) (M := A) (x := c - u * ϕ c) (y := b)).mpr
    (fun n => (sModEq_pow_iff I n _ _).mpr (key n))

/-- ★★★★★★**純抽象核 2(一意性)**。`I` 進 Hausdorff なだけで
半線型方程式 `c − u·ϕ(c) = b`(`u ∈ I`)の解は高々 1 つ。

★完備性は要らない。差 `d := c − c′` は `d = u·ϕ(d)` を満たすので
`d ∈ I^n` が全ての `n` で従い、Hausdorff で `d = 0`。 -/
theorem semilinear_solution_unique_of_isHausdorff {A : Type*} [CommRing A] (I : Ideal A)
    [IsHausdorff I A] (ϕ : A →+* A) (hϕ : ∀ x ∈ I, ϕ x ∈ I) {u : A} (hu : u ∈ I) {c c' : A}
    (h : c - u * ϕ c = c' - u * ϕ c') : c = c' := by
  have hd : c - c' = u * ϕ (c - c') := by
    rw [map_sub, mul_sub]
    linear_combination h
  have hmem : ∀ n : ℕ, c - c' ∈ I ^ n := by
    intro n
    induction n with
    | zero => simp
    | succ n ih =>
        rw [hd, pow_succ' I n]
        exact Ideal.mul_mem_mul hu (map_mem_pow_of_map_mem ϕ hϕ n ih)
  exact (IsHausdorff.eq_iff_smodEq (I := I) (M := A) (x := c) (y := c')).mpr
    (fun n => (sModEq_pow_iff I n c c').mpr (hmem n))

/-- ★★**既存の十分条件との整合(その 1・一般形)**。

`hsolve` の「解の集合」は Hausdorff な設定では高々 1 点である。
★すなわち `LubinTateEndoTwisted` が持つ 3 つの十分条件
(有限性 + 中山 / `ϕ = id` / 本ファイルの完備性)は、
**重なる範囲で必ず同じ `c` を返す**。 -/
theorem semilinear_solution_subsingleton {A : Type*} [CommRing A] (I : Ideal A) [IsHausdorff I A]
    (ϕ : A →+* A) (hϕ : ∀ x ∈ I, ϕ x ∈ I) {u : A} (hu : u ∈ I) (b : A) :
    Subsingleton {c : A // c - u * ϕ c = b} :=
  ⟨fun c c' => Subtype.ext
    (semilinear_solution_unique_of_isHausdorff I ϕ hϕ hu (c.2.trans c'.2.symm))⟩

/-- ★★**既存の十分条件との整合(その 2・`ϕ = id` の場合を明示)**。

`TwistedLT.ofUntwisted` は `ϕ = id` のとき `1 − u` が単数であることを使って
`c = (1−u)⁻¹b` と解く。Hausdorff ならば完備版の解 `c` はこれと一致する。 -/
theorem eq_of_one_sub_mul_eq_of_isHausdorff {A : Type*} [CommRing A] (I : Ideal A)
    [IsHausdorff I A] {u : A} (hu : u ∈ I) {b c w : A}
    (hc : c - u * (RingHom.id A) c = b) (hw : (1 - u) * w = b) : c = w := by
  refine semilinear_solution_unique_of_isHausdorff I (RingHom.id A) (fun _ hx => hx) hu ?_
  rw [hc]
  show b = w - u * w
  rw [← hw]
  ring

/-- ★**純抽象核 3**: 環同型は局所環の極大イデアルを保つ。

★`𝔪 = nonunits` と、mathlib の instance `isLocalHom_equiv`
(`Mathlib/Algebra/Group/Units/Equiv.lean`)だけを使う。
★これが `ϕ(𝔪) ⊆ 𝔪` の出どころである——付値も不分岐性も使っていない。 -/
theorem maximalIdeal_map_mem_of_ringEquiv {A : Type*} [CommRing A] [IsLocalRing A]
    (e : A ≃+* A) (x : A) (hx : x ∈ IsLocalRing.maximalIdeal A) :
    e x ∈ IsLocalRing.maximalIdeal A := by
  rw [IsLocalRing.mem_maximalIdeal, mem_nonunits_iff] at hx ⊢
  exact fun hu => hx (IsLocalHom.map_nonunit (f := e) x hu)

/-- 全単射な環準同型についての同じこと。 -/
theorem maximalIdeal_map_mem_of_bijective {A : Type*} [CommRing A] [IsLocalRing A]
    (ϕ : A →+* A) (hbij : Function.Bijective ϕ) (x : A)
    (hx : x ∈ IsLocalRing.maximalIdeal A) : ϕ x ∈ IsLocalRing.maximalIdeal A :=
  maximalIdeal_map_mem_of_ringEquiv (RingEquiv.ofBijective ϕ hbij) x hx

/-- 核 2 の、`ϕ` が環同型の場合。★`ϕ(𝔪) ⊆ 𝔪` が自動で埋まるうえ、
イデアルを `𝔪` に固定するので呼び出し側で `_` の推論が走らない
(`lean-idioms.md` #191: ここを `_` のままにすると `whnf` が heartbeat を焼く)。 -/
theorem semilinear_solution_unique_of_isHausdorff_ringEquiv {A : Type*} [CommRing A]
    [IsLocalRing A] [IsHausdorff (IsLocalRing.maximalIdeal A) A] (e : A ≃+* A)
    {u : A} (hu : u ∈ IsLocalRing.maximalIdeal A) {c c' : A}
    (h : c - u * e c = c' - u * e c') : c = c' :=
  semilinear_solution_unique_of_isHausdorff (IsLocalRing.maximalIdeal A) (e : A →+* A)
    (fun x hx => maximalIdeal_map_mem_of_ringEquiv e x hx) hu h

/-! ## 1. `hsolve` の完備版 -/

/-- ★★★★★★`TwistedLT` の `hsolve` の**第 3 の十分条件(完備性)**。

★`LubinTateEndoTwisted` にあった 2 つ(`hsolve_of_moduleFinite` = 有限性 + 中山、
`TwistedLT.ofUntwisted` = `ϕ = id`)はどちらも `L = K̂^ur` には当たらない
(`𝒪_{K̂^ur}` は `𝒪_K` 上有限でなく、`ϕ` は恒等でない)。★これが埋めた穴である。 -/
theorem hsolve_of_isAdicComplete {A : Type*} [CommRing A] [IsLocalRing A]
    [IsAdicComplete (IsLocalRing.maximalIdeal A) A] (ϕ : A →+* A)
    (hϕ : ∀ x ∈ IsLocalRing.maximalIdeal A, ϕ x ∈ IsLocalRing.maximalIdeal A) :
    ∀ u ∈ IsLocalRing.maximalIdeal A, ∀ b : A, ∃ c : A, c - u * ϕ c = b :=
  fun u hu b => exists_sub_mul_map_eq_of_isAdicComplete _ ϕ hϕ hu b

/-- `ϕ` が環同型なら `ϕ(𝔪) ⊆ 𝔪` は自動で埋まる。 -/
theorem hsolve_of_isAdicComplete_ringEquiv {A : Type*} [CommRing A] [IsLocalRing A]
    [IsAdicComplete (IsLocalRing.maximalIdeal A) A] (e : A ≃+* A) :
    ∀ u ∈ IsLocalRing.maximalIdeal A, ∀ b : A, ∃ c : A, c - u * e c = b :=
  hsolve_of_isAdicComplete (e : A →+* A) (fun x hx => maximalIdeal_map_mem_of_ringEquiv e x hx)

/-! ## 2. 標準 Lubin-Tate 級数(舞台が空でないことの材料)

原典の (3.2.1) を満たす最も簡単な級数 `f = πX + X^q`。 -/

/-- 標準 Lubin-Tate 級数 `f(X) = πX + X^q`。 -/
noncomputable def standardLubinTateSeries {A : Type*} [CommRing A] (π : A) (q : ℕ) :
    PowerSeries A := PowerSeries.C π * PowerSeries.X + PowerSeries.X ^ q

theorem constantCoeff_standardLubinTateSeries {A : Type*} [CommRing A] (π : A) {q : ℕ}
    (hq : q ≠ 0) : PowerSeries.constantCoeff (standardLubinTateSeries π q) = 0 := by
  simp [standardLubinTateSeries, hq]

/-- `f ≡ πX (mod deg 2)`。 -/
theorem coeff_one_standardLubinTateSeries {A : Type*} [CommRing A] (π : A) {q : ℕ}
    (hq : 2 ≤ q) : PowerSeries.coeff 1 (standardLubinTateSeries π q) = π := by
  rw [standardLubinTateSeries, map_add, PowerSeries.coeff_C_mul, PowerSeries.coeff_one_X,
    PowerSeries.coeff_X_pow]
  simp [show ¬ (1 = q) by omega]

/-- `f ≡ X^q (mod 𝔪)`。 -/
theorem map_residue_standardLubinTateSeries {A : Type*} [CommRing A] [IsLocalRing A] {π : A}
    (hπ : π ∈ IsLocalRing.maximalIdeal A) (q : ℕ) :
    PowerSeries.map (IsLocalRing.residue A) (standardLubinTateSeries π q) = PowerSeries.X ^ q := by
  rw [standardLubinTateSeries, map_add, map_mul, PowerSeries.map_C, PowerSeries.map_X, map_pow,
    PowerSeries.map_X, (IsLocalRing.residue_eq_zero_iff π).mpr hπ, map_zero, zero_mul, zero_add]

/-! ## 3. 舞台の構成子 -/

/-- ★★★★★**完備性から舞台を作る**。`hsolve` は `IsAdicComplete` が、
`ϕ(𝔪) ⊆ 𝔪` は `ϕ` が環同型であることが埋める。 -/
noncomputable def TwistedLT.ofAdicComplete {A : Type*} [CommRing A] [IsLocalRing A]
    [IsAdicComplete (IsLocalRing.maximalIdeal A) A] {pp ff : ℕ} (e : A ≃+* A)
    (hϕres : ∀ a : A, IsLocalRing.residue A (e a) = IsLocalRing.residue A a ^ pp ^ ff)
    {π π' : A} (hπmem : π ∈ IsLocalRing.maximalIdeal A)
    (hπ'max : IsLocalRing.maximalIdeal A = Ideal.span {π'})
    (f : PowerSeries A) (hf0 : PowerSeries.constantCoeff f = 0)
    (hf1 : PowerSeries.coeff 1 f = π)
    (hfres : PowerSeries.map (IsLocalRing.residue A) f = PowerSeries.X ^ pp ^ ff)
    (f' : PowerSeries A) (hf'0 : PowerSeries.constantCoeff f' = 0)
    (hf'1 : PowerSeries.coeff 1 f' = π')
    (hf'res : PowerSeries.map (IsLocalRing.residue A) f' = PowerSeries.X ^ pp ^ ff) :
    TwistedLT A pp ff where
  ϕ := (e : A →+* A)
  π := π
  π' := π'
  hπmem := hπmem
  hπ'max := hπ'max
  hϕres := hϕres
  hsolve := hsolve_of_isAdicComplete_ringEquiv e
  f := f
  hf0 := hf0
  hf1 := hf1
  hfres := hfres
  f' := f'
  hf'0 := hf'0
  hf'1 := hf'1
  hf'res := hf'res

/-! ## 4. 木の実物(`𝒪_{K̂^ur}` と算術 Frobenius)への当てはめ -/

variable {p : ℕ} [Fact p.Prime]

/-- ★★★★★★★★**`hsolve` の完備版が `K̂^ur` で実際に成り立つ**。

`σ` は任意の `unramGal K` でよい(剰余体での条件は要らない)。
★これが「完備版が空虚でない」ことの証拠である。 -/
theorem hsolve_unramifiedCompletionInt (K : PAdicLocalField p) (σ : unramGal K) :
    ∀ u ∈ IsLocalRing.maximalIdeal ↥(unramifiedCompletionInt K),
      ∀ b : ↥(unramifiedCompletionInt K),
        ∃ c : ↥(unramifiedCompletionInt K),
          c - u * unramGalCompletionInt K σ c = b := by
  haveI := isAdicComplete_unramifiedCompletionInt K
  exact hsolve_of_isAdicComplete_ringEquiv (unramGalCompletionInt K σ)

/-- 算術 Frobenius 版。 -/
theorem hsolve_arithFrobenius_unramifiedCompletionInt (K : PAdicLocalField p) :
    ∀ u ∈ IsLocalRing.maximalIdeal ↥(unramifiedCompletionInt K),
      ∀ b : ↥(unramifiedCompletionInt K),
        ∃ c : ↥(unramifiedCompletionInt K),
          c - u * unramGalCompletionInt K (arithFrobenius K) c = b :=
  hsolve_unramifiedCompletionInt K (arithFrobenius K)

/-- `K̂^ur` の解も一意である(`𝔪` 進完備なので Hausdorff)。 -/
theorem unique_hsolve_unramifiedCompletionInt (K : PAdicLocalField p) (σ : unramGal K)
    {u : ↥(unramifiedCompletionInt K)}
    (hu : u ∈ IsLocalRing.maximalIdeal ↥(unramifiedCompletionInt K))
    {c c' : ↥(unramifiedCompletionInt K)}
    (h : c - u * unramGalCompletionInt K σ c = c' - u * unramGalCompletionInt K σ c') :
    c = c' := by
  haveI := isAdicComplete_unramifiedCompletionInt K
  exact semilinear_solution_unique_of_isHausdorff_ringEquiv (unramGalCompletionInt K σ) hu h

/-- ★★★★★★**舞台を `A = 𝒪_{K̂^ur}` に当てはめる**。

`hϕres` は `residue_unramGalCompletionInt`(Λ6 の在庫)と
`residueCard_eq_pow`(`Nat.card 𝓀[K] = p^{f}`)から、
`hsolve` は完備版から、それぞれ埋まる。 -/
noncomputable def TwistedLT.ofUnramifiedCompletion (K : PAdicLocalField p) {σ : unramGal K}
    (hσ : ∀ w : ↥(unramifiedClosureInt K),
      unramGalResidue K σ (IsLocalRing.residue _ w)
        = (IsLocalRing.residue _ w) ^ (Nat.card (Valued.ResidueField K.carrier)))
    {π π' : ↥(unramifiedCompletionInt K)}
    (hπmem : π ∈ IsLocalRing.maximalIdeal ↥(unramifiedCompletionInt K))
    (hπ'max : IsLocalRing.maximalIdeal ↥(unramifiedCompletionInt K) = Ideal.span {π'})
    (f : PowerSeries ↥(unramifiedCompletionInt K))
    (hf0 : PowerSeries.constantCoeff f = 0)
    (hf1 : PowerSeries.coeff 1 f = π)
    (hfres : PowerSeries.map (IsLocalRing.residue _) f
      = PowerSeries.X ^ p ^ absoluteInertiaDegree K)
    (f' : PowerSeries ↥(unramifiedCompletionInt K))
    (hf'0 : PowerSeries.constantCoeff f' = 0)
    (hf'1 : PowerSeries.coeff 1 f' = π')
    (hf'res : PowerSeries.map (IsLocalRing.residue _) f'
      = PowerSeries.X ^ p ^ absoluteInertiaDegree K) :
    TwistedLT ↥(unramifiedCompletionInt K) p (absoluteInertiaDegree K) :=
  letI := isAdicComplete_unramifiedCompletionInt K
  TwistedLT.ofAdicComplete (unramGalCompletionInt K σ)
    (fun a => by
      rw [residue_unramGalCompletionInt K hσ a, residueCard_eq_pow K])
    hπmem hπ'max f hf0 hf1 hfres f' hf'0 hf'1 hf'res

/-- `2 ≤ q`(`q = Nat.card 𝓀[K]`)。剰余体は有限体なので位数 `≥ 2`。 -/
theorem two_le_residueCard_pow (K : PAdicLocalField p) : 2 ≤ p ^ absoluteInertiaDegree K := by
  rw [← residueCard_eq_pow]
  exact Finite.one_lt_card

/-- ★★★★★★★★★**舞台は空でない**。

`K̂^ur` の素元 `π′` を取り、`f = f′ = π′X + X^q`(`q = Nat.card 𝓀[K] = p^{f}`)と
算術 Frobenius で `TwistedLT ↥(unramifiedCompletionInt K) p (absoluteInertiaDegree K)` の
元が実際に作れる。

★これで `LubinTateEndoTwisted` の機械(`[θ]_{f,f′}`・関数等式・合成則)が
**`K̂^ur` の上で本当に動く**ことが確かめられた。 -/
theorem nonempty_twistedLT_unramifiedCompletionInt (K : PAdicLocalField p) :
    Nonempty (TwistedLT ↥(unramifiedCompletionInt K) p (absoluteInertiaDegree K)) := by
  obtain ⟨π, hπ0, hπmax⟩ := exists_uniformizer_carrierIntegers K
  set π' : ↥(unramifiedCompletionInt K) := uniformizerCompletionInt K π with hπ'def
  have hπ'max : IsLocalRing.maximalIdeal ↥(unramifiedCompletionInt K) = Ideal.span {π'} :=
    maximalIdeal_unramifiedCompletionInt_eq_span K hπ0 hπmax
  have hπ'mem : π' ∈ IsLocalRing.maximalIdeal ↥(unramifiedCompletionInt K) := by
    rw [hπ'max]
    exact Ideal.mem_span_singleton_self π'
  have hq2 : 2 ≤ p ^ absoluteInertiaDegree K := two_le_residueCard_pow K
  have hq0 : p ^ absoluteInertiaDegree K ≠ 0 := by omega
  refine ⟨TwistedLT.ofUnramifiedCompletion K (arithFrobenius_residue K) hπ'mem hπ'max
    (standardLubinTateSeries π' (p ^ absoluteInertiaDegree K))
    (constantCoeff_standardLubinTateSeries π' hq0)
    (coeff_one_standardLubinTateSeries π' hq2)
    (map_residue_standardLubinTateSeries hπ'mem _)
    (standardLubinTateSeries π' (p ^ absoluteInertiaDegree K))
    (constantCoeff_standardLubinTateSeries π' hq0)
    (coeff_one_standardLubinTateSeries π' hq2)
    (map_residue_standardLubinTateSeries hπ'mem _)⟩

end ABC3.Found.PGC
