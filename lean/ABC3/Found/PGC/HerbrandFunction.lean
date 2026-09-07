import ABC3.Found.PGC.FixedRingRamificationIndex

/-!
# Herbrand の定理と Herbrand 関数 φ_H(Yoshida 2008 Proposition 6.9)

典拠: T. Yoshida, *Local Class Field Theory via Lubin-Tate Theory*
(Ann. Fac. Sci. Toulouse 17-2, 2008; arXiv math/0606108) **Proposition 6.9**(物理 p.16)。
構造化済み原文は `ResearchPaper/1_Structured/Local Class Field Theory via Lubin-Tate Theory/`
`section-6.html` の `id="prop-6-9"`(`data-pdf-page="16"`,
`data-item="Proposition 6.9 (Herbrand)"`)。

原文 (Yoshida08 p.16):
> Proposition 6.9 (Herbrand). Define φ_H(n) := −1+1|H| _τ∈H min{i(τ), n+1} for n ∈ R[bb]_≥0. Also, for n ∈ R[bb]_≥0, define G_n := {σ ∈ G | i(σ) ≥ n + 1}, i.e. G_n = G_i if i ∈ Z[bb]_≥0 and n ∈ (i − 1, i]. Then G_nH/H = (G/H)_φ_H(n) for all n ∈ R[bb]_≥0.

★`pdftotext` は分数の横線と総和記号 `Σ` を出力しない(構造化 HTML の `data-txt` が
`1|H|` / 空文字で明示している)。読み下すと

  `φ_H(n) = −1 + (1/|H|) Σ_{τ∈H} min{i(τ), n+1}`,  `G_n = {σ ∈ G | i(σ) ≥ n+1}`。

設定は §6.2 (`#setup-6-2`):
> 6.2. The Hasse-Arf theorem. Let G = Gal(K′/K) with K′/K totally ramified as before,
> and let G ▷ H with G/H = Gal(K′′/K). For σ ∈ G, let σ[bar] = σH ∈ G/H be its image.

## 原典の Proof(`.txt` 1060–1065 行)と本ファイルの対応

> Proof. For σ ∈ G/H, replace σ by the element in σH which has the maximal value of i,
> and let i(σ) = m. Let τ ∈ H. If i(τ) ≥ m, then i(στ) ≥ m, hence i(στ) = m.
> If i(τ) < m, then i(τ) ≥ min{i(στ), i(σ−1)}, hence i(στ) = i(τ).
> Therefore i(στ) = min{i(τ), m}. Now the Lemma 6.8 gives i(σ[bar]) = φH(m − 1) + 1.
> Therefore, as φH is increasing, for n ∈ R≥0 we have
> σ ∈ GnH/H ⟺ m ≥ n + 1 ⟺ i(σ[bar]) ≥ φH(n) + 1 ⟺ σ ∈ (G/H)φH(n).

★★**この段落の第 1〜3 文には分岐も付値も Galois も出てこない。** 出てくるのは
「群 `G`、部分群 `H`、`i : G → Γ` が超距離不等式 `min{i(a), i(b)} ≤ i(ab)` と
`i(a⁻¹) = i(a)` を満たす」だけである。したがって §1 で**純群論**に切り出した
(`ultrametric_mul_eq_min` / `exists_max_on_coset` / `upperSubgroup`)。

## 本ファイルの構成

**§1 抽象核 —— 純群論**(分岐・付値・Galois の語彙が 1 つも出てこない)

* `ultrametric_le_of_mul_le` : 原文の `i(τ) ≥ min{i(στ), i(σ⁻¹)}` の帰結
  `i(στ) ≤ i(σ) → i(στ) ≤ i(τ)`。
* `ultrametric_mul_eq_min` : ★原文の `i(στ) = min{i(τ), m}` そのもの。
  ★仮定は `i(στ) ≤ i(σ)` の 1 本だけ —— 原文の 2 つの場合分け
  (`i(τ) ≥ m` / `i(τ) < m`)は不要だった。
* `exists_max_on_coset` : 原文の「replace σ by the element in σH which has the maximal
  value of i」。有限集合上で値が最大の元を取るだけ。
* `upperSubgroup` : 「上に閉じた述語 `P` について `{σ | P (i σ)}` は部分群」。
  ★`G_n` が部分群であることの根拠はこれ 1 本(超距離不等式 + `i(σ⁻¹) = i(σ)` + `P (i 1)`)。
* `mem_mul_subgroup_iff` : `σ ∈ S·H ↔ ∃ τ ∈ H, στ ∈ S`(商群を作らずに `G_nH/H` を書くため)。

**§2 抽象核 —— `ℕ∞` と `ℝ` の橋**

* `RealLeENat r x` : `r ≤ x`(`r : ℝ`, `x : ℕ∞`)。★`x = ⊤` は常に真。
* `truncENat x r` : `min{x, r}` の実数値版。★`x = ⊤` なら `r`。
* `realLeENat_natCast_iff` : ★★**`(m : ℝ) ≤ x ↔ (m : ℕ∞) ≤ x`** ——
  これが「実数添字」と「既存の ℕ 添字」を繋ぐ唯一の橋である。

**§3 `i` の超距離不等式**(具体層。★在庫に無かったので新たに立てた)

* `min_ramIndex_le_ramIndex_mul` : `min{i(σ), i(τ)} ≤ i(στ)`。
  段取りは `(στ)•α − α = σ•(τ•α − α) + (σ•α − α)` と `addVal_add` + `addVal_smul` だけ。
* `ramIndex_inv` : `i(σ⁻¹) = i(σ)`。

**§4 実数添字の下付き分岐群**

* `ramificationGroupReal α n` : ★原典の `G_n := {σ ∈ G | i(σ) ≥ n+1}` をそのまま定義にした。
* `ramificationGroupReal_eq_of_mem_Ioc` : ★★原典の "i.e. `G_n = G_i` if `i ∈ ℤ≥0` and
  `n ∈ (i − 1, i]`" そのもの。★**既存の ℕ 添字 `lowerRamificationGroup` と一致する。**
* `ramificationGroupReal_natCast` : その `n = m : ℕ` の場合。
  ★**`n < i(σ)`(既存、ℕ)と `n+1 ≤ i(σ)`(原典、ℝ)が同じであることを確かめた。**

**§5 Herbrand 関数 φ_H**

* `herbrandSum α H n := Σ_{τ∈H} min{i(τ), n+1}`(★**掛け算形。除算が無い**)。
* `herbrandPhi α H n := −1 + herbrandSum α H n / |H|`(★**除算はここだけ**)。
* `card_mul_herbrandPhi_add_one` : `|H| · (φ_H(n) + 1) = herbrandSum α H n`
  ——★これで下流は除算を見ずに済む。
* `strictMono_herbrandPhi` : ★原文の "as φ_H is increasing"。★**狭義**である
  (`τ = 1` の項が `min{⊤, n+1} = n+1` で狭義単調だから)。原文の最後の同値は
  `⟸` の向きで狭義性を使う。

**§6 Herbrand の定理**

* `herbrand_mem_iff` : `(∃ τ ∈ H, στ ∈ G_n) ↔ σ ∈ (G/H)_{φ_H(n)}`。
* `herbrand_coe_mul_coe_eq` : ★★★★**主結論** `G_n · H = (G/H)_{φ_H(n)} の引き戻し`。

## ★★商群 `G ⧸ H` を作っていない(Y9・Y10 と同じ)

原典は `σ̄ ∈ G/H` の `i(σ̄)` を `K′′` 上で測る。本ファイルは Y9 に倣い、
`G` を `C`(= `O_{K′′}` に対応する DVR)に直接作用させ、`ramIndex ϖ (σ : G)` と書く。
`H` が `C` に自明に作用する(`hHtriv`)ので、この関数は `G/H` を経由する。
したがって `G_nH/H = (G/H)_{φ_H(n)}` は、`G` の中の集合の等式

`(G_n : Set G) * (H : Set G) = ((G/H)_{φ_H(n)} の引き戻し : Set G)`

として述べてある。★**`MulSemiringAction (G ⧸ H) C` の構成は要らない。**
★同じ理由で `H` の正規性 `G ▷ H` も明示的には要らない。

## ★★除算をどこに置いたか(持ち場からの問い)

`φ_H` は本質的に有理数値なので、`lean-idioms.md` #102 の「除算を書かない」方針は
ここでは通らない。採った形は

* 値の型は **`ℝ`**(`ℚ` ではない)。理由: 原典が `n ∈ ℝ≥0` と実数の添字を要求しており、
  `φ_H` は `n` の関数として `ℝ → ℝ` でなければ `(G/H)_{φ_H(n)}` が書けない。
* **除算は `herbrandPhi` の定義式 1 箇所だけ**。証明本体はすべて掛け算形
  `card_mul_herbrandPhi_add_one` を経由する。Y9 の
  `card_mul_ramIndex_eq_sum_ramIndex`(掛け算形の Lemma 6.8)がそのまま噛み合った。
* `min{i(τ), n+1}` は `truncENat` で `ℝ` に落とす(★`ℕ∞` の引き算は使わない。#102)。
  `i(τ) = ⊤` のときだけ場合分けが要るが、`min` を取った後は常に有限なので安全である。

## ★★逸脱の記録

1. ★**原典は `n ∈ ℝ≥0` と非負に制限しているが、本ファイルの `ramificationGroupReal` と
   `herbrandPhi` は `n : ℝ` 全体で定義してある。** 制限しない方が定義も補題も短く、
   主結論 `herbrand_mem_iff` は `n < 0` でも真である(両辺とも `n` について反単調に
   広がるだけ)。★**弱めても強めてもいない。** `n ≥ 0` を落として偽になる主張は
   本ファイルには無い —— 唯一注意が要る `ramificationGroupReal_eq_of_mem_Ioc` は
   `m : ℕ` と `n ∈ ((m:ℝ)−1, m]` を**両方**仮定しており、そこで自動的に `n > −1` になる。
2. Lemma 6.8 は Y10 の `card_mul_ramIndex_eq_sum_ramIndex_of_fixedRing` を使う
   (`hram` を仮定に置かない形)。したがって本ファイルの仮定 `hcomp` `hHtriv` `hπ'`
   `hinj` `hfixC` `hres` `hAC` `hϖ` `hadj` `hfix` は**すべて Y10 のものをそのまま
   引き継いだ**ものであり、本ファイルが新たに追加した仮定は 1 つも無い。
   ★`hres`(`K′/K′′` の完全分岐)と `hAC`(`O ⊆ O_{K′′}`)は依然として仮定である
   (Y10 の報告どおり、上流がこれらを供給する小ノードを 1 つ必要とする)。
3. ★原典は `i(σ̄) = φ_H(m−1) + 1`(Lemma 6.8 の帰結)を**独立の等式として**書くが、
   本ファイルはそれを `herbrand_mem_iff` の証明の中の `hJval` として持っている。
   ★`m = ⊤`(すなわち `σ ∈ H`)のときこの等式は `⊤ = ⊤` でしか意味を持たないので、
   独立の宣言にすると `m ≠ ⊤` の仮定を外から与える必要が出る。主結論はその場合分けを
   内側で閉じている。

## 退化の自己検査

* ★**`i(1) = ⊤` を落とすと `φ_H` が壊れる**。`strictMono_herbrandSum` は
  `τ = 1` の項 `truncENat ⊤ (n+1) = n+1` が狭義単調であることだけに依っている
  (他の項は単調でしかない)。`ramIndex_one` を落とすと `φ_H` は狭義単調でなくなり、
  原文最後の同値の `⟸` の向きが出ない。
* ★**`H` が無限だと `Nat.card H = 0` で除算が壊れる**(`φ_H = −1` に潰れる)。
  `[Fintype ↥H]` を要求してある。Y9 の `pos_natCard_subgroup` と同じ穴。
* ★**「剰余類の中で `i` が最大の元」を取らずに任意の代表で議論すると偽**。
  `exists_max_on_coset` で取り替えている。取り替えを落とすと `i(στ) = min{i(τ), m}`
  (`ultrametric_mul_eq_min` の仮定 `i(στ) ≤ i(σ)`)が成り立たない。
* ★**`n ∈ ℝ` を `ℕ` に落とすと `(G/H)_{φ_H(n)}` が書けない** ——
  `φ_H(n)` は一般に整数でない(それが Herbrand 関数の要点である)。
* ★`ENat.mul_top` の穴: `h68` から `i_ϖ(σ) ≠ ⊤` を出すところで `0 < |H|` を使う。
  `|H| = 0` だと `|H| · i_ϖ(σ) = 0` で空虚に真になる。
-/

namespace ABC3.Found.PGC

open IsLocalRing IsDiscreteValuationRing
open scoped Pointwise

def herbrand_coe_mul_coe_eq.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 16, item := "Proposition 6.9 (Herbrand)",
    sectionId := "prop-6-9" }

/-! ## §1 抽象核 —— 純群論

★分岐・付値・Galois の語彙が 1 つも出てこない。原典の Proof の第 1〜3 文はここに尽きる。 -/

/-- **原文「If i(τ) < m, then i(τ) ≥ min{i(στ), i(σ⁻¹)}, hence i(στ) = i(τ)」の核** ——

超距離不等式 `min{f(a), f(b)} ≤ f(ab)` と `f(a⁻¹) = f(a)` だけから、
`f(στ) ≤ f(σ)` ならば `f(στ) ≤ f(τ)` が出る。

★段取り: `τ = σ⁻¹ · (στ)` に超距離不等式を当て、`min{f(σ⁻¹), f(στ)} = f(στ)` を使う。 -/
theorem ultrametric_le_of_mul_le {G Γ : Type*} [Group G] [LinearOrder Γ] {f : G → Γ}
    (hmul : ∀ a b : G, min (f a) (f b) ≤ f (a * b)) (hinv : ∀ a : G, f a⁻¹ = f a)
    {σ τ : G} (h : f (σ * τ) ≤ f σ) : f (σ * τ) ≤ f τ := by
  have h2 := hmul σ⁻¹ (σ * τ)
  rw [hinv, ← mul_assoc, inv_mul_cancel, one_mul, min_eq_right h] at h2
  exact h2

/-- ★★★**抽象核の本体 —— 原文「Therefore i(στ) = min{i(τ), m}」**。

`f : G → Γ` が超距離不等式と `f(a⁻¹) = f(a)` を満たし、`f(στ) ≤ f(σ)` なら
`f(στ) = min{f(τ), f(σ)}`。

★原文は `i(τ) ≥ m` と `i(τ) < m` の 2 つの場合に分けているが、
**場合分けは要らなかった** —— 必要なのは `f(στ) ≤ f(σ)`(剰余類の中で `σ` が最大)の 1 本だけ。 -/
theorem ultrametric_mul_eq_min {G Γ : Type*} [Group G] [LinearOrder Γ] {f : G → Γ}
    (hmul : ∀ a b : G, min (f a) (f b) ≤ f (a * b)) (hinv : ∀ a : G, f a⁻¹ = f a)
    {σ τ : G} (h : f (σ * τ) ≤ f σ) : f (σ * τ) = min (f τ) (f σ) :=
  le_antisymm (le_min (ultrametric_le_of_mul_le hmul hinv h) h)
    (by rw [min_comm]; exact hmul σ τ)

/-- **原文「replace σ by the element in σH which has the maximal value of i」** ——

有限部分群 `H` の剰余類 `σH` の中で `f` が最大になる代表 `σρ` が取れる。
★純粋に「有限集合の上で値が最大の元を取る」だけである。 -/
theorem exists_max_on_coset {G Γ : Type*} [Group G] [LinearOrder Γ] (f : G → Γ)
    (H : Subgroup G) [Fintype H] (σ : G) :
    ∃ ρ ∈ H, ∀ τ ∈ H, f (σ * τ) ≤ f (σ * ρ) := by
  obtain ⟨b, -, hb⟩ := Finset.exists_max_image (Finset.univ : Finset H)
    (fun τ : H => f (σ * (τ : G))) ⟨1, Finset.mem_univ _⟩
  exact ⟨b, b.2, fun τ hτ => hb ⟨τ, hτ⟩ (Finset.mem_univ _)⟩

/-- ★★**`G_n` が部分群であることの根拠**(抽象核) ——

`Γ` の上の**上に閉じた**述語 `P` と、超距離不等式 `min{f a, f b} ≤ f (ab)`・
`f(a⁻¹) = f(a)`・`P (f 1)` から、`{σ | P (f σ)}` は部分群になる。

★`mul_mem` は「`min (f a) (f b)` は `f a` か `f b` のどちらか」だけで出る
(`LinearOrder` であることが効く)。 -/
def upperSubgroup {G Γ : Type*} [Group G] [LinearOrder Γ] (f : G → Γ) (P : Γ → Prop)
    (hP : ∀ a b : Γ, a ≤ b → P a → P b)
    (hmul : ∀ a b : G, min (f a) (f b) ≤ f (a * b)) (hinv : ∀ a : G, f a⁻¹ = f a)
    (hone : P (f 1)) : Subgroup G where
  carrier := {σ | P (f σ)}
  mul_mem' {a b} ha hb := by
    refine hP _ _ (hmul a b) ?_
    rcases min_cases (f a) (f b) with ⟨h, _⟩ | ⟨h, _⟩ <;> rw [h] <;> assumption
  one_mem' := hone
  inv_mem' {a} ha := by rw [Set.mem_setOf_eq, hinv]; exact ha

@[simp] theorem mem_upperSubgroup {G Γ : Type*} [Group G] [LinearOrder Γ] {f : G → Γ}
    {P : Γ → Prop} {hP hmul hinv hone} {σ : G} :
    σ ∈ upperSubgroup f P hP hmul hinv hone ↔ P (f σ) := Iff.rfl

/-- **`G_nH/H` を商群を作らずに書くための言い換え** —— `σ ∈ S·H ↔ ∃ τ ∈ H, στ ∈ S`。 -/
theorem mem_mul_subgroup_iff {G : Type*} [Group G] (S H : Subgroup G) (σ : G) :
    σ ∈ (S : Set G) * (H : Set G) ↔ ∃ τ ∈ H, σ * τ ∈ S := by
  rw [Set.mem_mul]
  constructor
  · rintro ⟨x, hx, y, hy, rfl⟩
    exact ⟨y⁻¹, H.inv_mem hy, by simpa using hx⟩
  · rintro ⟨τ, hτ, h⟩
    exact ⟨σ * τ, h, τ⁻¹, H.inv_mem hτ, by group⟩

/-! ## §2 抽象核 —— `ℕ∞` と `ℝ` の橋

★`i` は `ℕ∞` に値を取る(`ramIndex`)が、原典の添字 `n` は実数である。
両者を繋ぐ語彙をここで 2 つだけ用意する。★`ℕ∞` の側では引き算も除算もしない(#102)。 -/

/-- `r ≤ x`(`r : ℝ`, `x : ℕ∞`)。★`x = ⊤`(すなわち `i(1)`)のときは常に真。 -/
def RealLeENat (r : ℝ) (x : ℕ∞) : Prop := x = ⊤ ∨ r ≤ (x.toNat : ℝ)

theorem RealLeENat.top {r : ℝ} : RealLeENat r ⊤ := Or.inl rfl

@[simp] theorem realLeENat_coe {r : ℝ} {k : ℕ} : RealLeENat r (k : ℕ∞) ↔ r ≤ (k : ℝ) := by
  simp [RealLeENat]

/-- `r` について反変・`x` について共変。★`G_n` が部分群であること・`n` について
反単調であることの両方がこの 1 本から出る。 -/
theorem RealLeENat.mono {r s : ℝ} {x y : ℕ∞} (hrs : s ≤ r) (hxy : x ≤ y)
    (h : RealLeENat r x) : RealLeENat s y := by
  rcases eq_or_ne y ⊤ with rfl | hy
  · exact Or.inl rfl
  have hx : x ≠ ⊤ := fun hc => hy (top_le_iff.1 (hc ▸ hxy))
  rcases h with h | h
  · exact absurd h hx
  refine Or.inr (hrs.trans (h.trans ?_))
  exact_mod_cast Nat.cast_le.2 (by exact_mod_cast (ENat.coe_toNat hx) ▸ (ENat.coe_toNat hy) ▸ hxy)

/-- ★★**実数添字と ℕ 添字を繋ぐ唯一の橋** —— `(m : ℝ) ≤ x ↔ (m : ℕ∞) ≤ x`。 -/
theorem realLeENat_natCast_iff {m : ℕ} {x : ℕ∞} : RealLeENat (m : ℝ) x ↔ (m : ℕ∞) ≤ x := by
  rcases eq_or_ne x ⊤ with rfl | hx
  · simp [RealLeENat.top]
  · rw [← ENat.coe_toNat hx, realLeENat_coe]
    exact_mod_cast Iff.rfl

/-- 原文の `min{i(τ), n+1}` の実数値版。★`x = ⊤` のときは `r`(切り詰めが効く)。
★`min` を取った後は常に有限なので `⊤` を安全に消せる。 -/
noncomputable def truncENat (x : ℕ∞) (r : ℝ) : ℝ := if x = ⊤ then r else min (x.toNat : ℝ) r

@[simp] theorem truncENat_top (r : ℝ) : truncENat ⊤ r = r := by simp [truncENat]

@[simp] theorem truncENat_coe (k : ℕ) (r : ℝ) : truncENat (k : ℕ∞) r = min (k : ℝ) r := by
  simp [truncENat]

/-- ★`φ_H` の単調性はこれを各項に当てるだけで出る。 -/
theorem monotone_truncENat (x : ℕ∞) : Monotone (truncENat x) := by
  intro r s hrs
  unfold truncENat; split
  · exact hrs
  · exact min_le_min le_rfl hrs

/-- `ℕ∞` で取った `min{x, M}` を `ℕ` 経由で `ℝ` に落とすと `truncENat` に一致する。
★Lemma 6.8(`ℕ∞` の等式)を `φ_H`(`ℝ` の等式)に翻訳する要。 -/
theorem toNat_min_natCast (x : ℕ∞) (M : ℕ) :
    (((min x (M : ℕ∞)).toNat : ℕ) : ℝ) = truncENat x (M : ℝ) := by
  cases x with
  | top => simp
  | coe k =>
      simp only [truncENat_coe]
      rcases le_total k M with h | h
      · rw [min_eq_left (by exact_mod_cast h : (k : ℕ∞) ≤ (M : ℕ∞)), ENat.toNat_coe,
          min_eq_left (by exact_mod_cast h : (k : ℝ) ≤ (M : ℝ))]
      · rw [min_eq_right (by exact_mod_cast h : (M : ℕ∞) ≤ (k : ℕ∞)), ENat.toNat_coe,
          min_eq_right (by exact_mod_cast h : (M : ℝ) ≤ (k : ℝ))]

theorem min_natCast_ne_top (x : ℕ∞) (M : ℕ) : min x (M : ℕ∞) ≠ ⊤ :=
  ne_top_of_le_ne_top (ENat.coe_ne_top M) (min_le_right _ _)

/-! ## §3 `i` の超距離不等式(具体層)

★在庫に無かったので新たに立てた。原典の Proof が `i(τ) ≥ min{i(στ), i(σ⁻¹)}` の形で使う。 -/

/-- ★★**`i` の超距離不等式** `min{i(σ), i(τ)} ≤ i(στ)`。

段取りは環の恒等式 `(στ)•α − α = σ•(τ•α − α) + (σ•α − α)` に
`addVal` の劣加法性(`addVal_add`)と群作用不変性(`addVal_smul`)を当てるだけ。
★`α` が素元であることは要らない。 -/
theorem min_ramIndex_le_ramIndex_mul {B : Type*} [CommRing B] [IsDomain B]
    [IsDiscreteValuationRing B] {G : Type*} [Group G] [MulSemiringAction G B]
    (α : B) (σ τ : G) :
    min (ramIndex α σ) (ramIndex α τ) ≤ ramIndex α (σ * τ) := by
  have key : (σ * τ) • α - α = σ • (τ • α - α) + (σ • α - α) := by
    rw [mul_smul, smul_sub]; ring
  calc min (ramIndex α σ) (ramIndex α τ)
      = min (addVal B (σ • (τ • α - α))) (addVal B (σ • α - α)) := by
        simp only [ramIndex]; rw [addVal_smul]; exact min_comm _ _
    _ ≤ addVal B (σ • (τ • α - α) + (σ • α - α)) := addVal_add
    _ = ramIndex α (σ * τ) := by rw [ramIndex, key]

/-- **`i(σ⁻¹) = i(σ)`** —— `σ⁻¹•α − α = −σ⁻¹•(σ•α − α)`。 -/
theorem ramIndex_inv {B : Type*} [CommRing B] [IsDomain B]
    [IsDiscreteValuationRing B] {G : Type*} [Group G] [MulSemiringAction G B]
    (α : B) (σ : G) : ramIndex α σ⁻¹ = ramIndex α σ := by
  have key : σ⁻¹ • α - α = -(σ⁻¹ • (σ • α - α)) := by
    rw [smul_sub, ← mul_smul, inv_mul_cancel, one_smul]; ring
  rw [ramIndex, ramIndex, key, AddValuation.map_neg, addVal_smul]

/-! ## §4 実数添字の下付き分岐群 -/

/-- ★★**原典の `G_n := {σ ∈ G | i(σ) ≥ n + 1}`(`n` は実数)**。

部分群であることは §1 の `upperSubgroup` に、超距離不等式(§3)・`i(σ⁻¹) = i(σ)`(§3)・
`i(1) = ⊤`(`ramIndex_one`)を渡すだけで出る。

★逸脱: 原典は `n ∈ ℝ≥0` だが、ここは `n : ℝ` 全体で定義してある(冒頭「逸脱の記録 (1)」)。 -/
noncomputable def ramificationGroupReal {B : Type*} [CommRing B] [IsDomain B]
    [IsDiscreteValuationRing B] {G : Type*} [Group G] [MulSemiringAction G B]
    (α : B) (n : ℝ) : Subgroup G :=
  upperSubgroup (ramIndex α) (RealLeENat (n + 1))
    (fun _ _ hab h => h.mono le_rfl hab)
    (min_ramIndex_le_ramIndex_mul α) (ramIndex_inv α)
    (by rw [ramIndex_one]; exact RealLeENat.top)

@[simp] theorem mem_ramificationGroupReal {B : Type*} [CommRing B] [IsDomain B]
    [IsDiscreteValuationRing B] {G : Type*} [Group G] [MulSemiringAction G B]
    {α : B} {n : ℝ} {σ : G} :
    σ ∈ ramificationGroupReal α n ↔ RealLeENat (n + 1) (ramIndex α σ) := Iff.rfl

/-- ★★★**原典「i.e. `G_n = G_i` if `i ∈ ℤ≥0` and `n ∈ (i − 1, i]`」** ——

実数添字の `G_n` は、`n ∈ ((m:ℝ)−1, m]` なる `m : ℕ` に対して
**既存の ℕ 添字 `lowerRamificationGroup B G m` と一致する**。

★★これが「実数へ延ばした定義が下流とずれていない」ことの検算である。
既存の在庫は `σ ∈ G_m ↔ (m : ℕ∞) < i(σ)` の形(`ℕ` の狭義不等号)で、
原典は `i(σ) ≥ n+1` の形(`ℝ` の非狭義不等号)。両者が一致することをここで確かめた。 -/
theorem ramificationGroupReal_eq_of_mem_Ioc {A B : Type*} [CommRing A] [CommRing B] [Algebra A B]
    [IsDomain B] [IsDiscreteValuationRing B] {G : Type*} [Group G] [MulSemiringAction G B]
    [SMulCommClass G A B] {α : B} (huni : maximalIdeal B = Ideal.span {α})
    (hadj : Algebra.adjoin A ({α} : Set B) = ⊤) {n : ℝ} (m : ℕ)
    (h1 : (m : ℝ) - 1 < n) (h2 : n ≤ (m : ℝ)) :
    ramificationGroupReal α n = lowerRamificationGroup B G m := by
  ext σ
  rw [mem_ramificationGroupReal, mem_lowerRamificationGroup_iff_lt_ramIndex (A := A) huni hadj]
  rcases eq_or_ne (ramIndex α σ) ⊤ with hx | hx
  · rw [hx]
    exact iff_of_true RealLeENat.top (Ne.lt_top (ENat.coe_ne_top m))
  · rw [← ENat.coe_toNat hx, realLeENat_coe]
    constructor
    · intro h
      have hlt : m < (ramIndex α σ).toNat := by
        by_contra hc
        have : ((ramIndex α σ).toNat : ℝ) ≤ (m : ℝ) := by exact_mod_cast Nat.not_lt.1 hc
        linarith
      exact_mod_cast hlt
    · intro h
      have hlt : m < (ramIndex α σ).toNat := by exact_mod_cast h
      have : (m : ℝ) + 1 ≤ ((ramIndex α σ).toNat : ℝ) := by exact_mod_cast hlt
      linarith

/-- ★**整数添字での一致**(`ramificationGroupReal_eq_of_mem_Ioc` の `n = m` の場合)。 -/
theorem ramificationGroupReal_natCast {A B : Type*} [CommRing A] [CommRing B] [Algebra A B]
    [IsDomain B] [IsDiscreteValuationRing B] {G : Type*} [Group G] [MulSemiringAction G B]
    [SMulCommClass G A B] {α : B} (huni : maximalIdeal B = Ideal.span {α})
    (hadj : Algebra.adjoin A ({α} : Set B) = ⊤) (m : ℕ) :
    ramificationGroupReal α (m : ℝ) = lowerRamificationGroup B G m :=
  ramificationGroupReal_eq_of_mem_Ioc huni hadj m (by linarith) le_rfl

/-- `G_n` は `n` について反単調(既存の `lowerRamificationGroup_antitone` の実数版)。 -/
theorem ramificationGroupReal_antitone {B : Type*} [CommRing B] [IsDomain B]
    [IsDiscreteValuationRing B] {G : Type*} [Group G] [MulSemiringAction G B] (α : B) :
    Antitone (ramificationGroupReal (G := G) α) := by
  intro a b hab σ hσ
  exact RealLeENat.mono (by linarith) le_rfl hσ

/-! ## §5 Herbrand 関数 `φ_H`

★★ここで初めて除算が避けられない(`φ_H` は本質的に有理数値)。
除算は `herbrandPhi` の定義式 1 箇所に閉じ込め、証明本体は掛け算形
`card_mul_herbrandPhi_add_one` を経由する。 -/

/-- 原文 `Σ_{τ∈H} min{i(τ), n+1}` —— ★**掛け算形。除算が無い**。
`|H| · (φ_H(n) + 1)` に等しい(`card_mul_herbrandPhi_add_one`)。 -/
noncomputable def herbrandSum {B : Type*} [CommRing B] [IsDomain B] [IsDiscreteValuationRing B]
    {G : Type*} [Group G] [MulSemiringAction G B] (α : B) (H : Subgroup G) [Fintype H]
    (n : ℝ) : ℝ :=
  ∑ τ : H, truncENat (ramIndex α (τ : G)) (n + 1)

/-- ★★**原典の Herbrand 関数** `φ_H(n) := −1 + (1/|H|) Σ_{τ∈H} min{i(τ), n+1}`。

★★**本ファイルで除算が現れるのはここだけである。** -/
noncomputable def herbrandPhi {B : Type*} [CommRing B] [IsDomain B] [IsDiscreteValuationRing B]
    {G : Type*} [Group G] [MulSemiringAction G B] (α : B) (H : Subgroup G) [Fintype H]
    (n : ℝ) : ℝ :=
  -1 + herbrandSum α H n / (Nat.card H : ℝ)

/-- ★★**掛け算形** `|H| · (φ_H(n) + 1) = Σ_{τ∈H} min{i(τ), n+1}`。
下流(Herbrand 本体)はこれだけを使い、除算を見ない。 -/
theorem card_mul_herbrandPhi_add_one {B : Type*} [CommRing B] [IsDomain B]
    [IsDiscreteValuationRing B] {G : Type*} [Group G] [MulSemiringAction G B] (α : B)
    (H : Subgroup G) [Fintype H] (n : ℝ) :
    (Nat.card H : ℝ) * (herbrandPhi α H n + 1) = herbrandSum α H n := by
  have h : (Nat.card H : ℝ) ≠ 0 := Nat.cast_ne_zero.2 (Nat.card_pos (α := H)).ne'
  simp only [herbrandPhi]
  field_simp
  ring

/-- ★**原文「as φ_H is increasing」の核** —— `Σ_{τ∈H} min{i(τ), n+1}` は `n` について
**狭義**単調増加。★狭義性は `τ = 1` の項 `min{⊤, n+1} = n+1` だけから来る
(`ramIndex_one` を落とすと崩れる)。 -/
theorem strictMono_herbrandSum {B : Type*} [CommRing B] [IsDomain B]
    [IsDiscreteValuationRing B] {G : Type*} [Group G] [MulSemiringAction G B] (α : B)
    (H : Subgroup G) [Fintype H] : StrictMono (herbrandSum α H) := by
  intro a b hab
  refine Finset.sum_lt_sum (fun τ _ => monotone_truncENat _ (by linarith))
    ⟨1, Finset.mem_univ _, ?_⟩
  rw [Subgroup.coe_one, ramIndex_one, truncENat_top, truncENat_top]
  linarith

/-- ★★**原文「as φ_H is increasing」** —— `φ_H` は狭義単調増加。
★原文最後の同値の `⟸` の向きは**狭義**性を使う(単調だけでは出ない)。 -/
theorem strictMono_herbrandPhi {B : Type*} [CommRing B] [IsDomain B]
    [IsDiscreteValuationRing B] {G : Type*} [Group G] [MulSemiringAction G B] (α : B)
    (H : Subgroup G) [Fintype H] : StrictMono (herbrandPhi α H) := by
  intro a b hab
  have hpos : (0 : ℝ) < (Nat.card H : ℝ) := by exact_mod_cast Nat.card_pos (α := H)
  simp only [herbrandPhi, add_lt_add_iff_left]
  exact (div_lt_div_iff_of_pos_right hpos).mpr (strictMono_herbrandSum α H hab)

/-! ## §6 Herbrand の定理(Yoshida 2008 Proposition 6.9) -/

/-- ★★★★**Yoshida 2008 Proposition 6.9 (Herbrand)** —— 所属の形。

`σ ∈ G_nH/H ⟺ σ̄ ∈ (G/H)_{φ_H(n)}`、すなわち

`(∃ τ ∈ H, στ ∈ G_n) ↔ i_ϖ(σ) ≥ φ_H(n) + 1`。

原文 (Yoshida08 p.16):
> Then G_nH/H = (G/H)_φ_H(n) for all n ∈ R[bb]_≥0.

* `C` は `O_{K′′}` に対応する DVR。`H` は `C` に自明に作用する(`hHtriv`)ので
  `ramIndex ϖ (·)` は `G/H` を経由する。★**商群 `G ⧸ H` は作らない。**
* 仮定はすべて Y10 の `card_mul_ramIndex_eq_sum_ramIndex_of_fixedRing`(仮定なしの
  Lemma 6.8)のものをそのまま引き継いだ。★**本定理が新たに足した仮定は 1 つも無い。**

★段取り(原文どおり):
1. `exists_max_on_coset` で剰余類 `σH` の中の `i` 最大の代表 `σρ` を取る(`m := i(σρ)`)。
2. `ultrametric_mul_eq_min` で `i(σρτ) = min{i(τ), m}`。
3. Lemma 6.8 で `|H| · i_ϖ(σ) = Σ_{τ∈H} min{i(τ), m}`。
4. `m` が有限なら右辺は `herbrandSum π' H (m−1)`、すなわち `i_ϖ(σ) = φ_H(m−1)+1`。
5. `strictMono_herbrandPhi` で `φ_H(n) ≤ φ_H(m−1) ⟺ n ≤ m−1 ⟺ m ≥ n+1`。
6. `m = ⊤` の場合(`σ ∈ H`)は両辺とも真(`ramIndex_one`)。 -/
theorem herbrand_mem_iff {A B : Type*} [CommRing A] [CommRing B] [Algebra A B] [IsDomain B]
    [IsDiscreteValuationRing B] {G : Type*} [Group G] [MulSemiringAction G B]
    [SMulCommClass G A B] [FaithfulSMul G B] {H : Subgroup G} [Fintype H] {π' π'' : B}
    {C : Type*} [CommRing C] [IsDomain C] [IsDiscreteValuationRing C] [Algebra C B]
    [MulSemiringAction G C]
    (hcomp : ∀ (ρ : G) (c : C), algebraMap C B (ρ • c) = ρ • algebraMap C B c)
    (hHtriv : ∀ ρ ∈ H, ∀ c : C, ρ • c = c)
    (hπ' : Irreducible π')
    (hinj : Function.Injective (algebraMap C B))
    (hfixC : ∀ b : B, (∀ ρ ∈ H, ρ • b = b) → ∃ c : C, algebraMap C B c = b)
    (hres : ∀ b : B, ∃ c : C, b - algebraMap C B c ∈ maximalIdeal B)
    (hAC : ∀ a : A, ∃ c : C, algebraMap C B c = algebraMap A B a)
    {ϖ : C} (hϖ : algebraMap C B ϖ = π'')
    (hadj : Algebra.adjoin A ({π'} : Set B) = ⊤)
    (hfix : ∀ y : B, (∀ ρ ∈ H, ρ • y = y) → y ∈ Algebra.adjoin A ({π''} : Set B))
    (n : ℝ) (σ : G) :
    (∃ τ ∈ H, σ * τ ∈ ramificationGroupReal π' n)
      ↔ σ ∈ ramificationGroupReal ϖ (herbrandPhi π' H n) := by
  classical
  have hcardpos : 0 < Nat.card H := Nat.card_pos
  have hcardR : (0 : ℝ) < (Nat.card H : ℝ) := by exact_mod_cast hcardpos
  obtain ⟨ρ, hρ, hmax⟩ := exists_max_on_coset (fun s : G => ramIndex π' s) H σ
  have hcoset : ∀ τ : H, ramIndex π' (σ * ρ * (τ : G))
      = min (ramIndex π' (τ : G)) (ramIndex π' (σ * ρ)) := by
    intro τ
    refine ultrametric_mul_eq_min (min_ramIndex_le_ramIndex_mul π') (ramIndex_inv π') ?_
    rw [mul_assoc]
    exact hmax _ (H.mul_mem hρ τ.2)
  have hLHS : (∃ τ ∈ H, σ * τ ∈ ramificationGroupReal π' n)
      ↔ RealLeENat (n + 1) (ramIndex π' (σ * ρ)) := by
    refine ⟨?_, fun h => ⟨ρ, hρ, h⟩⟩
    rintro ⟨τ, hτ, hmem⟩
    exact RealLeENat.mono le_rfl (hmax τ hτ) hmem
  have hϖeq : ramIndex ϖ (σ * ρ) = ramIndex ϖ σ := by
    show addVal C ((σ * ρ) • ϖ - ϖ) = addVal C (σ • ϖ - ϖ)
    rw [mul_smul, hHtriv ρ hρ ϖ]
  have h68 : (Nat.card H : ℕ∞) * ramIndex ϖ σ
      = ∑ τ : H, min (ramIndex π' (τ : G)) (ramIndex π' (σ * ρ)) := by
    rw [← hϖeq, card_mul_ramIndex_eq_sum_ramIndex_of_fixedRing (A := A) hcomp hHtriv hπ' hinj
      hfixC hres hAC hϖ hadj hfix (σ * ρ)]
    exact Finset.sum_congr rfl fun τ _ => hcoset τ
  rw [hLHS, mem_ramificationGroupReal]
  rcases eq_or_ne (ramIndex π' (σ * ρ)) ⊤ with hm | hm
  · -- `m = ⊤`、すなわち `σρ = 1`(`σ ∈ H`)。両辺とも真。
    have h1 : σ * ρ = 1 :=
      eq_one_of_smul_eq_of_adjoin_eq_top hadj ((ramIndex_eq_top_iff π' (σ * ρ)).1 hm)
    have h2 : ramIndex ϖ σ = ⊤ := by rw [← hϖeq, h1, ramIndex_one]
    rw [hm, h2]
    exact iff_of_true RealLeENat.top RealLeENat.top
  · -- `m = M : ℕ` が有限の場合。
    have hMeq : ramIndex π' (σ * ρ) = ((ramIndex π' (σ * ρ)).toNat : ℕ∞) :=
      (ENat.coe_toNat hm).symm
    set M : ℕ := (ramIndex π' (σ * ρ)).toNat with hMdef
    have hsum : ∑ τ : H, min (ramIndex π' (τ : G)) (ramIndex π' (σ * ρ))
        = ((∑ τ : H, (min (ramIndex π' (τ : G)) (M : ℕ∞)).toNat : ℕ) : ℕ∞) := by
      rw [hMeq, Nat.cast_sum]
      exact Finset.sum_congr rfl fun τ _ => (ENat.coe_toNat (min_natCast_ne_top _ _)).symm
    rw [hsum] at h68
    have hJne : ramIndex ϖ σ ≠ ⊤ := by
      intro hc
      rw [hc, ENat.mul_top (by exact_mod_cast hcardpos.ne')] at h68
      exact (ENat.coe_ne_top _) h68.symm
    have hJeq : ramIndex ϖ σ = ((ramIndex ϖ σ).toNat : ℕ∞) := (ENat.coe_toNat hJne).symm
    set J : ℕ := (ramIndex ϖ σ).toNat with hJdef
    have hnat : Nat.card H * J = ∑ τ : H, (min (ramIndex π' (τ : G)) (M : ℕ∞)).toNat := by
      rw [hJeq] at h68
      exact_mod_cast h68
    have hreal : (Nat.card H : ℝ) * (J : ℝ) = herbrandSum π' H ((M : ℝ) - 1) := by
      have h1 : ((Nat.card H * J : ℕ) : ℝ)
          = ((∑ τ : H, (min (ramIndex π' (τ : G)) (M : ℕ∞)).toNat : ℕ) : ℝ) := by
        exact_mod_cast hnat
      rw [Nat.cast_mul, Nat.cast_sum] at h1
      rw [h1]
      unfold herbrandSum
      refine Finset.sum_congr rfl fun τ _ => ?_
      rw [show (M : ℝ) - 1 + 1 = (M : ℝ) by ring]
      exact toNat_min_natCast _ _
    -- ★原文「the Lemma 6.8 gives i(σ[bar]) = φH(m − 1) + 1」
    have hJval : (J : ℝ) = herbrandPhi π' H ((M : ℝ) - 1) + 1 := by
      rw [← card_mul_herbrandPhi_add_one π' H ((M : ℝ) - 1)] at hreal
      exact mul_left_cancel₀ hcardR.ne' hreal
    rw [hMeq, hJeq, realLeENat_coe, realLeENat_coe, hJval]
    constructor
    · intro h
      have hle := (strictMono_herbrandPhi π' H).monotone (show n ≤ (M : ℝ) - 1 by linarith)
      linarith
    · intro h
      have hle : herbrandPhi π' H n ≤ herbrandPhi π' H ((M : ℝ) - 1) := by linarith
      have := (strictMono_herbrandPhi π' H).le_iff_le.1 hle
      linarith

/-- ★★★★**Yoshida 2008 Proposition 6.9 (Herbrand)** —— 集合の等式の形。

原文 (Yoshida08 p.16):
> Then G_nH/H = (G/H)_φ_H(n) for all n ∈ R[bb]_≥0.

`G ⧸ H` を作らずに、`G` の中の集合の等式

`(G_n : Set G) * (H : Set G) = ((G/H)_{φ_H(n)} の引き戻し : Set G)`

として述べた。右辺は `H` が `C` に自明に作用する(`hHtriv`)ので `H` 剰余類の合併であり、
商 `G/H` の分岐群の引き戻しに他ならない。 -/
theorem herbrand_coe_mul_coe_eq {A B : Type*} [CommRing A] [CommRing B] [Algebra A B]
    [IsDomain B] [IsDiscreteValuationRing B] {G : Type*} [Group G] [MulSemiringAction G B]
    [SMulCommClass G A B] [FaithfulSMul G B] {H : Subgroup G} [Fintype H] {π' π'' : B}
    {C : Type*} [CommRing C] [IsDomain C] [IsDiscreteValuationRing C] [Algebra C B]
    [MulSemiringAction G C]
    (hcomp : ∀ (ρ : G) (c : C), algebraMap C B (ρ • c) = ρ • algebraMap C B c)
    (hHtriv : ∀ ρ ∈ H, ∀ c : C, ρ • c = c)
    (hπ' : Irreducible π')
    (hinj : Function.Injective (algebraMap C B))
    (hfixC : ∀ b : B, (∀ ρ ∈ H, ρ • b = b) → ∃ c : C, algebraMap C B c = b)
    (hres : ∀ b : B, ∃ c : C, b - algebraMap C B c ∈ maximalIdeal B)
    (hAC : ∀ a : A, ∃ c : C, algebraMap C B c = algebraMap A B a)
    {ϖ : C} (hϖ : algebraMap C B ϖ = π'')
    (hadj : Algebra.adjoin A ({π'} : Set B) = ⊤)
    (hfix : ∀ y : B, (∀ ρ ∈ H, ρ • y = y) → y ∈ Algebra.adjoin A ({π''} : Set B))
    (n : ℝ) :
    ((ramificationGroupReal (G := G) π' n : Subgroup G) : Set G) * (H : Set G)
      = ((ramificationGroupReal (G := G) ϖ (herbrandPhi π' H n) : Subgroup G) : Set G) := by
  ext σ
  rw [mem_mul_subgroup_iff, SetLike.mem_coe]
  exact herbrand_mem_iff (A := A) hcomp hHtriv hπ' hinj hfixC hres hAC hϖ hadj hfix n σ

end ABC3.Found.PGC
