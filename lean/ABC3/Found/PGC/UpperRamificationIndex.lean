import ABC3.Found.PGC.TameQuotientTower
import ABC3.Found.PGC.UpperRamificationGroup

/-!
# 上付き分岐群の指数 `|G/G^m| ∣ (q − 1)q^{m−1}`(Yoshida 2008 Corollary 6.13 (iii))

★★★**本ノードは Corollary 6.13 の (iii) のみ**を担当する。
**(i)(ii) は `Found/PGC/UpperRamificationGroup.lean`**(Y18)にある
(`upperRamification_coe_mul_coe_eq` / `upperRamificationGroup_eq_bot_of_two_quotients`)。
本ファイルは**既存の `Found/PGC/*.lean` を 1 行も書き換えない**(import のみ)。

典拠: T. Yoshida, *Local Class Field Theory via Lubin-Tate Theory*
(Ann. Fac. Sci. Toulouse 17-2, 2008; arXiv math/0606108) **Corollary 6.13**(物理 p.17)。
構造化済み原文は
`ResearchPaper/1_Structured/Local Class Field Theory via Lubin-Tate Theory/section-6.html` の
`id="cor-6-13"`(`data-pdf-page="17"`, `data-item="Corollary 6.13"`)。

原文 (Yoshida08 p.17):
> Corollary 6.13. (i) If G ⊲ H, then G^mH/H = (G/H)^m for all m ∈ R[bb]_≥0. (ii) Let K′/K and K′′/K be two Galois extensions with K′K′′/K totally ramified. If Gal(K′/K)^m = Gal(K′′/K)^m = {id} for m ∈ R[bb]_≥0, then Gal(K′K′′/K)^m = {id}. (iii) Let G be abelian. Then |G/G^m| divides (q − 1)q^m−1 for m ∈ Z[bb]_≥0.

★逐語の `(q − 1)q^m−1` は `pdftotext` が上付きの `m−1` を落とした形で、
原典は `(q − 1)q^{m−1}` である。

原典の Proof (iii)(`.txt` 1140–1143 行):
> (iii): If n − 1 < φ^{−1}_G(m) ≤ n for n ∈ Z≥0, then G^m = G_n. Consider G_i for integers
> 1 ≤ i ≤ n. Then, by Theorem 6.11, G_{i−1} ≠ G_i can only happen when φ_G(i − 1) ∈ Z,
> and as 0 ≤ φ_G(i − 1) ≤ φ_G(n − 1) < m, at most m − 1 times for i > 1.
> By Proposition 6.2, |G_{i−1}/G_i| divides q − 1 when i = 1 and q when i > 1.

## ★★これが今できるようになった理由

Y18 は (iii) を落としていた。理由は「Theorem 6.11 の**一般の可換 `G`** に対する形が要る」で、
その時点の `HasseArf.lean` は巡回群・馴分岐商までだったからである。
2026-09-07 に `TameQuotientTower.lean`(Y23)が
`exists_natCast_herbrandPhiGroup_of_abelian_of_setup` を**仮定ゼロ**で閉じたので解禁された。
★**この Hasse-Arf はそのまま使えた**(仮定の付け足しも読み替えも要らなかった)。

## ★★設計 —— まず抽象核(§1)、それから具体層(§2–§3)

原文の証明を読むと、**分岐・付値・Galois が 1 語も出てこない部分**が 3 つに割れる。

| § | 宣言 | 語彙 |
|---|---|---|
| §1 | `index_eq_prod_relIndex` | ★**純群論**。減少列の望遠鏡積 `[G:F n] = ∏ [F i : F (i+1)]` |
| §1 | `prod_dvd_pow_card_of_eq_one_outside` | ★**純算術**。`Finset` 上の積と `∣` だけ |
| §1 | `card_le_of_injOn_mem_Icc` | ★**純組合せ**。`{1,…,M}` への単射から `card ≤ M` |
| §1 | `index_dvd_mul_pow_of_jumps` | ★**純群論**。上の 3 本の合流 |
| §2 | `relIndex_lowerRamificationGroup_one_dvd` / `..._succ_dvd` | Prop 6.2 の 2 本の単射を `Nat.card` に落とすだけ |
| §3 | `index_upperRamificationGroup_dvd` ほか | 代入するだけ |

★`index_dvd_mul_pow_of_jumps` には **`B` も付値も出てこない**。
「望遠鏡積で位数を評価し、跳びの個数で `q` の冪を抑える」は群の減少列の一般論である。
★★**具体層は §1 に `F i := G_i`、`g i := ⌊φ_G(i)⌋₊` を代入するだけ**になった。

## ★本体の見当は当たった

段取りの見当(`θ_0` で `q−1` / `θ_n` で `q` / 跳びは高々 `m−1` 回 / 望遠鏡積)は
**4 つとも当たった**。差分は 2 つだけである。

1. ★**`|G_{i−1}/G_i|` は `Subgroup.relIndex` で書くのが安い。**
   `Subgroup.relIndex H K = Nat.card (K ⧸ H.subgroupOf K)` は `Subgroup.index_eq_card` で
   出るので、Y3 の `thetaMulQuot` / `thetaAddQuot` の定義域とそのまま合う。
   望遠鏡積は mathlib の `Subgroup.relIndex_mul_index` 1 本で、帰納法 4 行で終わる。
2. ★**跳びの数え上げに `φ_G` の「値」は要らない。** 要るのは
   「跳ぶ `i` に対して `⌊φ_G(i)⌋₊` が `{1,…,m−1}` に入り、`i ↦ ⌊φ_G(i)⌋₊` が単射」だけで、
   単射性は `φ_G` の狭義単調性(Y18 `strictMono_herbrandPhiGroup`)から出る。
   原文の `0 ≤ φ_G(i−1) ≤ φ_G(n−1) < m` は、この 2 つに分けると
   **`φ_G(0) = 0` と `φ_G(ψ_G(m)) = m` の 2 回の狭義単調性**に帰着する。

## ★★mathlib を先に引いた

★Y21 が「無い」と報告した 4 本を Y22 が mathlib で見つけた直後なので、
**新しい `def` / `instance` を書く前に `.cache/mathlib-index.txt` を名前空間で grep した**。
引き当ては次のとおりで、★**本ファイルが新たに定義する `def` / `instance` は 1 つも無い**。

| 要るもの | mathlib の宣言 | 引き当てた grep |
|---|---|---|
| `H.relIndex K * K.index = H.index` | `Subgroup.relIndex_mul_index` | `grep -n "relIndex_mul_index" .cache/mathlib-index.txt` |
| `H.relIndex H = 1` | `Subgroup.relIndex_self` | 同上(`relIndex` で 1 回) |
| `H.index = Nat.card (G ⧸ H)` | `Subgroup.index_eq_card` | 同上 |
| 単射準同型で `Nat.card` が割る | `Subgroup.card_dvd_of_injective` | `grep -n "card_dvd_of_injective"` |
| `Nat.card αˣ = Nat.card α − 1` | `Nat.card_units` | `grep -n "card_units"` |
| `∏ f ∣ ∏ g` | `Finset.prod_dvd_prod_of_dvd` | `grep -n "prod_dvd_prod_of_dvd"` |
| `InjOn` から `card ≤ card` | `Finset.card_le_card_of_injOn` | `grep -n "card_le_card_of_injOn"` |

★**`Nat.card_units` は `[GroupWithZero α]` だけで有限性を要求しない。**
おかげで剰余体の有限性を仮定に足さずに済んだ(下の「逸脱の記録 5」)。

## 逸脱の記録

1. **`m : ℕ` で述べ、`q^{m−1}` の `−1` は ℕ の切り詰め引き算**にした。
   原文は `m ∈ ℤ≥0` だが、`m = 0` では `(q−1)q^{−1}` が整数でないので字面のままでは
   主張が意味を持たない。ℕ の切り詰めでは `q^{0−1} = q^0 = 1` なので右辺は `q − 1` になり、
   左辺は `|G/G^0| = 1`(`index_upperRamificationGroup_zero`)なので**主張は真のまま**である。
   ★`m ≥ 1` の場合は字面どおり(`natCard_quot_upperRamificationGroup_succ_dvd`)。
2. **`q := Nat.card (ResidueField B)`**、すなわち `K′` の剰余体の位数とした。原文の `q` は
   `K` の剰余体の位数だが、仮定 `hresA`(`K′/K` が完全分岐 ⇒ 剰余体が伸びない)により
   両者は一致する。★`q` は `|G|` ではない。
3. **`G` 可換は `∀ x y : G, x * y = y * x` の命題**で渡す(先行ノード Y22・Y23 の流儀。
   `CommGroup` 構造だと具体層に当たらない)。
4. **底の設定は Y23 の Hasse-Arf(仮定ゼロ版)と同じ**である
   (`hπ'` / `hresA` / `hadj` / `hAinj` / `h0` と剰余標数 `p`)。本ノードが**足した仮定は無い**。
5. **剰余体の有限性を仮定していない。** 無限なら `Nat.card (ResidueField B) = 0` で
   右辺が `0` になり主張は自明に真。原典の設定(局所体)では有限で、そのとき非自明である。
6. **`upperRamificationGroup` は `m : ℝ` 全体で定義されている**(Y18 の逸脱 1)。本ファイルは
   `m : ℕ` を `((m : ℕ) : ℝ)` で流し込むだけで、`m ≥ 0` の制限は自動的に満たされる。

## 退化の自己検査

* ★★**`G` 可換を落とすと偽**。可換性は Theorem 6.11(Hasse-Arf)の仮定であり、
  これを落とすと「跳びは `φ_G` が整数値をとる点でしか起きない」が崩れ、
  跳びの回数を `m − 1` で抑える段が消える(`q` の冪の指数が抑えられない)。
* ★**`h0 : G_0 = ⊤`(完全分岐)を落とすと望遠鏡積の先頭が合わない。**
  `index_eq_prod_relIndex` は `F 0 = ⊤` を使って `[G : F n] = ∏_{i<n} [F i : F (i+1)]` を出す。
  `G_0 ≠ G` だと左辺に `|G/G_0|` の分だけ余分が出る。
* ★**`m` の `≥ 0`**: `m : ℕ` なので自動。`ψ_G(m) ≥ 0`(`herbrandPsiGroup_nonneg`)は
  `n := ⌈ψ_G(m)⌉₊` が `(n : ℝ) − 1 < ψ_G(m)` を満たすために使う。
* ★**`m = 0` の退化**: 逸脱 1 のとおり。`ψ_G(0) = 0` なので `n = 0`、
  `G^0 = G_0 = ⊤` で `|G/G^0| = 1`。積は空積 `1` になり、`q^{0−1} = 1` と辻褄が合う。
* ★**`i = 0` の項だけ `q − 1`、`i ≥ 1` の項が `q`** という非対称は Prop 6.2 の
  2 本の単射(`θ_0 : G_0/G_1 ↪ 𝓀ˣ` と `θ_i : G_i/G_{i+1} ↪ 𝓀⁺`)の行き先の違いそのもので、
  §1 の抽象核では `e`(先頭)と `q`(残り)という 2 つの独立なパラメータに分けてある。
* ★`ℕ∞` の切り詰め引き算・除算は**使っていない**(`lean-idioms.md` #102)。
-/

namespace ABC3.Found.PGC

open IsLocalRing IsDiscreteValuationRing

/-! ## §1 抽象核 —— 純群論・純算術

★**分岐・付値・Galois が 1 語も出てこない。** `Subgroup` と `Finset` と `ℕ` だけ。 -/

/-- ★★**抽象核(純群論)** —— 減少列の**望遠鏡積**
`[G : F n] = ∏_{i<n} [F i : F (i+1)]`。

`F 0 = ⊤` と `F (i+1) ≤ F i` だけを使う。mathlib の
`Subgroup.relIndex_mul_index : H ≤ K → H.relIndex K * K.index = H.index` の帰納。
★正規性も有限性も要らない(`Subgroup.index` は無限のとき `0` になる規約で、
そのときも等式は成り立つ)。 -/
theorem index_eq_prod_relIndex {G : Type*} [Group G] (F : ℕ → Subgroup G)
    (hle : ∀ i, F (i + 1) ≤ F i) (h0 : F 0 = ⊤) :
    ∀ n : ℕ, (F n).index = ∏ i ∈ Finset.range n, (F (i + 1)).relIndex (F i)
  | 0 => by simp [h0]
  | (n + 1) => by
      rw [Finset.prod_range_succ, ← index_eq_prod_relIndex F hle h0 n, mul_comm,
        Subgroup.relIndex_mul_index (hle n)]

/-- ★**抽象核(純算術)** —— `t` の外では `f i = 1` なら `∏_{i∈s} f i ∣ q ^ |t|`。

★原文の「`i > 1` では高々 `m − 1` 回」を支える形。`s ∩ t` に落として
`Finset.prod_dvd_prod_of_dvd` を 1 回叩くだけで、`q ≠ 0` も要らない。 -/
theorem prod_dvd_pow_card_of_eq_one_outside {ι : Type*} [DecidableEq ι] {s t : Finset ι}
    {f : ι → ℕ} {q : ℕ} (hdvd : ∀ i ∈ s, f i ∣ q) (h1 : ∀ i ∈ s, i ∉ t → f i = 1) :
    ∏ i ∈ s, f i ∣ q ^ t.card := by
  have heq : ∏ i ∈ s ∩ t, f i = ∏ i ∈ s, f i :=
    Finset.prod_subset Finset.inter_subset_left
      (fun x hx hnx => h1 x hx (fun hxt => hnx (Finset.mem_inter.2 ⟨hx, hxt⟩)))
  rw [← heq]
  refine dvd_trans (Finset.prod_dvd_prod_of_dvd _ _
    (fun i hi => hdvd i (Finset.mem_of_mem_inter_left hi))) ?_
  rw [Finset.prod_const]
  exact pow_dvd_pow q (Finset.card_le_card Finset.inter_subset_right)

/-- ★**抽象核(純組合せ)** —— `{1,…,M}` への単射があれば `|t| ≤ M`。

★原文の `0 ≤ φ_G(i−1) ≤ φ_G(n−1) < m` から「高々 `m − 1` 回」を出す段はこれだけ。 -/
theorem card_le_of_injOn_mem_Icc {ι : Type*} {t : Finset ι} {g : ι → ℕ} {M : ℕ}
    (hmem : ∀ i ∈ t, 1 ≤ g i ∧ g i ≤ M) (hinj : Set.InjOn g t) : t.card ≤ M := by
  have h := Finset.card_le_card_of_injOn (t := Finset.Icc 1 M) g
    (fun i hi => Finset.mem_Icc.2 (hmem i hi)) hinj
  simpa using h

/-- ★★★**抽象核(純群論)—— Corollary 6.13 (iii) の骨**。

減少列 `F : ℕ → Subgroup G`(`F 0 = ⊤`)に対して、

* 先頭の相対指数 `[F 0 : F 1]` が `e` を割り、
* `i ≥ 1` の相対指数 `[F i : F (i+1)]` が `q` を割り、
* `1 ≤ i < n` で `F i ≠ F (i+1)` となる `i` が「`{1,…,M}` に単射に入る」

なら `[G : F n] ∣ e * q ^ M`。

★★**分岐も付値も Galois も出てこない。** 具体層は `F i := G_i`、`e := q − 1`、
`M := m − 1`、`g i := ⌊φ_G(i)⌋₊` を代入するだけである。 -/
theorem index_dvd_mul_pow_of_jumps {G : Type*} [Group G] (F : ℕ → Subgroup G)
    (hle : ∀ i, F (i + 1) ≤ F i) (h0 : F 0 = ⊤) {e q M n : ℕ} {g : ℕ → ℕ}
    (he : (F 1).relIndex (F 0) ∣ e)
    (hq : ∀ i, 1 ≤ i → (F (i + 1)).relIndex (F i) ∣ q)
    (hbd : ∀ i, 1 ≤ i → i < n → F i ≠ F (i + 1) → 1 ≤ g i ∧ g i ≤ M)
    (hinj : ∀ i j, 1 ≤ i → i < n → F i ≠ F (i + 1) → 1 ≤ j → j < n → F j ≠ F (j + 1) →
      g i = g j → i = j) :
    (F n).index ∣ e * q ^ M := by
  classical
  rw [index_eq_prod_relIndex F hle h0 n]
  match n with
  | 0 => simp
  | (n' + 1) =>
    rw [Finset.prod_range_succ']
    set t : Finset ℕ := (Finset.range n').filter (fun i => F (i + 1) ≠ F (i + 2)) with ht
    have hprod : ∏ i ∈ Finset.range n', (F (i + 1 + 1)).relIndex (F (i + 1)) ∣ q ^ t.card := by
      refine prod_dvd_pow_card_of_eq_one_outside (fun i _ => hq (i + 1) (by omega)) ?_
      intro i hi hnit
      have hEq : F (i + 1) = F (i + 2) := by
        by_contra hc
        exact hnit (Finset.mem_filter.2 ⟨hi, hc⟩)
      rw [show i + 1 + 1 = i + 2 from rfl, ← hEq]
      exact Subgroup.relIndex_self _
    have hcard : t.card ≤ M := by
      refine card_le_of_injOn_mem_Icc (g := fun i => g (i + 1)) ?_ ?_
      · intro i hi
        rw [ht, Finset.mem_filter, Finset.mem_range] at hi
        exact hbd (i + 1) (by omega) (by omega) hi.2
      · intro i hi j hj hij
        simp only [ht, Finset.coe_filter, Set.mem_setOf_eq, Finset.mem_range] at hi hj
        have := hinj (i + 1) (j + 1) (by omega) (by omega) hi.2 (by omega) (by omega) hj.2 hij
        omega
    exact mul_dvd_mul (hprod.trans (pow_dvd_pow q hcard)) he |>.trans (by rw [mul_comm])

/-! ## §2 具体層 —— Proposition 6.2 の 2 本の単射を `relIndex` に落とす

★Y3 の `thetaMulQuot_injective` / `thetaAddQuot_injective` に
`Subgroup.card_dvd_of_injective` を当てるだけ。
★`Subgroup.relIndex H K = Nat.card (K ⧸ H.subgroupOf K)` は
`Subgroup.index_eq_card` で出て、`theta*Quot` の定義域とそのまま一致する。 -/

section Theta

variable {A B : Type*} [CommRing A] [CommRing B] [Algebra A B] [IsDomain B]
  [IsLocalRing B] {G : Type*} [Group G] [MulSemiringAction G B] [SMulCommClass G A B]

/-- ★★**Prop 6.2 前半の帰結** —— `[G_0 : G_1]` は `q − 1` を割る。

原文 (Yoshida08 p.17) の Corollary 6.13 (iii) の証明:
> By Proposition 6.2, |G_{i−1}/G_i| divides q − 1 when i = 1

★`Nat.card_units : Nat.card αˣ = Nat.card α − 1` は `[GroupWithZero α]` だけで
成り立つので、剰余体の**有限性を仮定しなくてよい**。 -/
theorem relIndex_lowerRamificationGroup_one_dvd
    {α : B} (hα : maximalIdeal B = Ideal.span {α}) (hα0 : α ≠ 0)
    (hadj : Algebra.adjoin A ({α} : Set B) = ⊤) :
    (lowerRamificationGroup B G 1).relIndex (lowerRamificationGroup B G 0)
      ∣ Nat.card (ResidueField B) - 1 := by
  rw [← Nat.card_units (α := ResidueField B),
    show (lowerRamificationGroup B G 1).relIndex (lowerRamificationGroup B G 0)
      = Nat.card (lowerRamificationGroup B G 0 ⧸
          (lowerRamificationGroup B G 1).subgroupOf (lowerRamificationGroup B G 0)) from
      Subgroup.index_eq_card _]
  exact Subgroup.card_dvd_of_injective _ (thetaMulQuot_injective (A := A) hα hα0 hadj)

/-- ★★**Prop 6.2 後半の帰結** —— `i ≥ 1` なら `[G_i : G_{i+1}]` は `q` を割る。

原文 (Yoshida08 p.17) の Corollary 6.13 (iii) の証明:
> and q when i > 1

★行き先は `Multiplicative (ResidueField B)`(剰余体の加法群を乗法的に見たもの)で、
`Nat.card (Multiplicative X) = Nat.card X` は型シノニムなので `show` 1 行で移せる。 -/
theorem relIndex_lowerRamificationGroup_succ_dvd
    {α : B} (hα : maximalIdeal B = Ideal.span {α}) (hα0 : α ≠ 0)
    (hadj : Algebra.adjoin A ({α} : Set B) = ⊤) {i : ℕ} (hi : 1 ≤ i) :
    (lowerRamificationGroup B G (i + 1)).relIndex (lowerRamificationGroup B G i)
      ∣ Nat.card (ResidueField B) := by
  rw [show (lowerRamificationGroup B G (i + 1)).relIndex (lowerRamificationGroup B G i)
      = Nat.card (lowerRamificationGroup B G i ⧸
          (lowerRamificationGroup B G (i + 1)).subgroupOf (lowerRamificationGroup B G i)) from
      Subgroup.index_eq_card _]
  show _ ∣ Nat.card (Multiplicative (ResidueField B))
  exact Subgroup.card_dvd_of_injective _ (thetaAddQuot_injective (A := A) hα hα0 hadj hi)

end Theta

/-! ## §3 Corollary 6.13 (iii)

★底の設定は Y23 の `exists_natCast_herbrandPhiGroup_of_abelian_of_setup` と**同じ**である
(本ノードが足した仮定は無い)。 -/

section Main

variable {A B : Type*} [CommRing A] [IsDomain A] [IsDiscreteValuationRing A] [CommRing B]
  [Algebra A B] [IsDomain B] [IsDiscreteValuationRing B] [IsNoetherian A B]
  {G : Type*} [Group G] [MulSemiringAction G B] [Fintype G] [SMulCommClass G A B]
  [FaithfulSMul G B] [Fintype ↥(lowerRamificationGroup B G 1)]

omit [FaithfulSMul G B] [Fintype ↥(lowerRamificationGroup B G 1)] in
/-- ★**`m = 0` の退化** —— `|G/G^0| = 1`。

原文は `m ∈ Z≥0` で `(q−1)q^{m−1}` と書くが、`m = 0` では `(q−1)q^{−1}` は整数ではない。
★`G^0 = G_0 = ⊤`(Y18 `upperRamificationGroup_zero`)なので左辺は `1` で、
どの右辺の読み方でも主張は真である。 -/
theorem index_upperRamificationGroup_zero {α : B}
    (hα : maximalIdeal B = Ideal.span {α}) :
    (upperRamificationGroup G α (0 : ℝ)).index = 1 := by
  rw [upperRamificationGroup_zero hα]
  exact Subgroup.index_top

/-- ★★**Hasse-Arf の消費** —— 跳ぶ `i` では `⌊φ_G(i)⌋₊` が `φ_G(i)` そのもの。

原文 (Yoshida08 p.17):
> by Theorem 6.11, G_i−1 = G_i can only happen when φ_G(i − 1) ∈ Z

★逐語の `G_i−1 = G_i` は `pdftotext` が `≠` の斜線を落とした形である
(原典は `G_{i−1} ≠ G_i`)。★添字の波括弧も `pdftotext` には出ない
——[[pdftotext-drops-negation]]。

★Y23 の `exists_natCast_herbrandPhiGroup_of_abelian_of_setup`(仮定ゼロの Hasse-Arf)を
**そのまま**使う。`∃ k : ℕ, φ_G(i) = k` から `⌊φ_G(i)⌋₊ = k` を読むだけで、
以降の数え上げでは `k` を名指しせず `⌊·⌋₊` で通せる。 -/
theorem natCast_floor_herbrandPhiGroup_of_jump
    (p : ℕ) [Fact p.Prime] [CharP (ResidueField B) p]
    {π' : B} (hπ' : Irreducible π')
    (hresA : ∀ b : B, ∃ a : A, b - algebraMap A B a ∈ maximalIdeal B)
    (hadj : Algebra.adjoin A ({π'} : Set B) = ⊤)
    (hAinj : Function.Injective (algebraMap A B))
    (h0 : lowerRamificationGroup B G 0 = ⊤) (habel : ∀ x y : G, x * y = y * x)
    {i : ℕ} (hne : lowerRamificationGroup B G i ≠ lowerRamificationGroup B G (i + 1)) :
    ((⌊herbrandPhiGroup G π' (i : ℝ)⌋₊ : ℕ) : ℝ) = herbrandPhiGroup G π' (i : ℝ) := by
  obtain ⟨k, hk⟩ := exists_natCast_herbrandPhiGroup_of_abelian_of_setup (A := A) p hπ' hresA hadj
    hAinj h0 habel hne
  rw [hk, Nat.floor_natCast]

def index_upperRamificationGroup_dvd.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 17, item := "Corollary 6.13", sectionId := "cor-6-13" }

/-- ★★★★★**Yoshida 2008 Corollary 6.13 (iii)** ——
`G` が可換なら `|G/G^m|` は `(q − 1)q^{m−1}` を割る。

原文 (Yoshida08 p.17):
> Corollary 6.13. (i) If G ⊲ H, then G^mH/H = (G/H)^m for all m ∈ R[bb]_≥0. (ii) Let K′/K and K′′/K be two Galois extensions with K′K′′/K totally ramified. If Gal(K′/K)^m = Gal(K′′/K)^m = {id} for m ∈ R[bb]_≥0, then Gal(K′K′′/K)^m = {id}. (iii) Let G be abelian. Then |G/G^m| divides (q − 1)q^m−1 for m ∈ Z[bb]_≥0.

★★**本宣言は (iii) のみ**。(i)(ii) は `Found/PGC/UpperRamificationGroup.lean`。

`q := Nat.card (ResidueField B)`(逸脱 2)、`q^{m−1}` の `−1` は ℕ の切り詰め引き算
(逸脱 1。`m = 0` では右辺が `q − 1` になり、左辺は `1` なので真)。

★★**仮定は Y23 の Hasse-Arf(仮定ゼロ版)と同じ**で、本ノードが足したものは無い:

* `hπ'` : `π′` は `𝒪_{K′}` の素元。
* `hresA` : `K′/K` が完全分岐(剰余体が伸びない)。★原文の `q` が `K′` 側でも `q` である根拠。
* `hadj` : `𝒪_{K′} = 𝒪_K[π′]`(原典 Lemma 5.11)。
* `hAinj` : `𝒪_K → 𝒪_{K′}` が単射。
* `h0` : `G_0 = G`(`K′/K` が完全分岐)。★望遠鏡積の先頭が合うために要る。
* `habel` : 原文の `Let G be abelian`。★これを落とすと Hasse-Arf が使えず偽。
* `[CharP (ResidueField B) p]` : 剰余標数 `p`。 -/
theorem index_upperRamificationGroup_dvd
    (p : ℕ) [Fact p.Prime] [CharP (ResidueField B) p]
    {π' : B} (hπ' : Irreducible π')
    (hresA : ∀ b : B, ∃ a : A, b - algebraMap A B a ∈ maximalIdeal B)
    (hadj : Algebra.adjoin A ({π'} : Set B) = ⊤)
    (hAinj : Function.Injective (algebraMap A B))
    (h0 : lowerRamificationGroup B G 0 = ⊤) (habel : ∀ x y : G, x * y = y * x)
    (m : ℕ) :
    (upperRamificationGroup G π' (m : ℝ)).index
      ∣ (Nat.card (ResidueField B) - 1) * Nat.card (ResidueField B) ^ (m - 1) := by
  have hα : maximalIdeal B = Ideal.span {π'} := (irreducible_iff_uniformizer π').1 hπ'
  have hα0 : π' ≠ 0 := hπ'.ne_zero
  -- 原文「If n − 1 < φ^{−1}_G(m) ≤ n for n ∈ Z≥0, then G^m = G_n」を `n := ⌈ψ_G(m)⌉₊` で使う。
  have hψnn : 0 ≤ herbrandPsiGroup G π' (m : ℝ) :=
    herbrandPsiGroup_nonneg hα (Nat.cast_nonneg m)
  have h2 : herbrandPsiGroup G π' (m : ℝ) ≤ (⌈herbrandPsiGroup G π' (m : ℝ)⌉₊ : ℝ) :=
    Nat.le_ceil _
  have h1 : ((⌈herbrandPsiGroup G π' (m : ℝ)⌉₊ : ℕ) : ℝ) - 1 < herbrandPsiGroup G π' (m : ℝ) := by
    have h := Nat.ceil_lt_add_one hψnn
    linarith
  rw [upperRamificationGroup_eq_lowerRamificationGroup (A := A) hα hadj _ h1 h2]
  refine index_dvd_mul_pow_of_jumps (fun i => lowerRamificationGroup B G i)
    (fun i => lowerRamificationGroup_antitone B G (Nat.le_succ i)) h0
    (g := fun i => ⌊herbrandPhiGroup G π' (i : ℝ)⌋₊)
    (relIndex_lowerRamificationGroup_one_dvd (A := A) hα hα0 hadj)
    (fun i hi => relIndex_lowerRamificationGroup_succ_dvd (A := A) hα hα0 hadj hi) ?_ ?_
  · -- 原文「0 ≤ φ_G(i−1) ≤ φ_G(n−1) < m ⇒ 高々 m − 1 回」
    intro i hi1 hin hne
    have hfl := natCast_floor_herbrandPhiGroup_of_jump (A := A) p hπ' hresA hadj hAinj h0 habel hne
    have hpos : (0 : ℝ) < herbrandPhiGroup G π' (i : ℝ) := by
      have hlt := strictMono_herbrandPhiGroup (G := G) π'
        (show (0 : ℝ) < (i : ℝ) by exact_mod_cast hi1)
      rwa [herbrandPhiGroup_zero hα] at hlt
    have hltm : herbrandPhiGroup G π' (i : ℝ) < (m : ℝ) := by
      have hiψ : (i : ℝ) < herbrandPsiGroup G π' (m : ℝ) := by
        have hle : ((i : ℕ) : ℝ) + 1 ≤ ((⌈herbrandPsiGroup G π' (m : ℝ)⌉₊ : ℕ) : ℝ) := by
          exact_mod_cast hin
        linarith
      have hlt := strictMono_herbrandPhiGroup (G := G) π' hiψ
      rwa [herbrandPhiGroup_herbrandPsiGroup hα] at hlt
    refine ⟨?_, ?_⟩
    · rw [Nat.one_le_iff_ne_zero]
      intro hc
      rw [hc, Nat.cast_zero] at hfl
      linarith
    · have hcast : ((⌊herbrandPhiGroup G π' (i : ℝ)⌋₊ : ℕ) : ℝ) < ((m : ℕ) : ℝ) := by
        rw [hfl]; exact hltm
      have : ⌊herbrandPhiGroup G π' (i : ℝ)⌋₊ < m := by exact_mod_cast hcast
      omega
  · -- 跳ぶ `i` たちは `φ_G` の狭義単調性で分離される。
    intro i j hi1 hin hnei hj1 hjn hnej hgij
    have hfi := natCast_floor_herbrandPhiGroup_of_jump (A := A) p hπ' hresA hadj hAinj h0 habel hnei
    have hfj := natCast_floor_herbrandPhiGroup_of_jump (A := A) p hπ' hresA hadj hAinj h0 habel hnej
    have hphi : herbrandPhiGroup G π' (i : ℝ) = herbrandPhiGroup G π' (j : ℝ) := by
      rw [← hfi, ← hfj, hgij]
    have := (strictMono_herbrandPhiGroup (G := G) π').injective hphi
    exact_mod_cast this

def natCard_quot_upperRamificationGroup_dvd.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 17, item := "Corollary 6.13", sectionId := "cor-6-13" }

/-- ★★★★**Corollary 6.13 (iii) の `|G/G^m|` の形**(原文の記法どおり)。

原文 (Yoshida08 p.17):
> (iii) Let G be abelian. Then |G/G^m| divides (q − 1)q^m−1 for m ∈ Z[bb]_≥0.

★上の `index` 版を `Subgroup.index_eq_card` で言い換えただけ。 -/
theorem natCard_quot_upperRamificationGroup_dvd
    (p : ℕ) [Fact p.Prime] [CharP (ResidueField B) p]
    {π' : B} (hπ' : Irreducible π')
    (hresA : ∀ b : B, ∃ a : A, b - algebraMap A B a ∈ maximalIdeal B)
    (hadj : Algebra.adjoin A ({π'} : Set B) = ⊤)
    (hAinj : Function.Injective (algebraMap A B))
    (h0 : lowerRamificationGroup B G 0 = ⊤) (habel : ∀ x y : G, x * y = y * x)
    (m : ℕ) :
    Nat.card (G ⧸ upperRamificationGroup G π' (m : ℝ))
      ∣ (Nat.card (ResidueField B) - 1) * Nat.card (ResidueField B) ^ (m - 1) := by
  rw [← Subgroup.index_eq_card]
  exact index_upperRamificationGroup_dvd (A := A) p hπ' hresA hadj hAinj h0 habel m

/-- ★★**`m ≥ 1` の形**(切り詰め引き算が出てこない版)。

`m = k + 1` と書くと右辺は `(q − 1) * q ^ k` で、原文の `(q − 1)q^{m−1}` の
**字面どおり**になる。★逸脱 1 が効くのは `m = 0` のときだけである。 -/
theorem natCard_quot_upperRamificationGroup_succ_dvd
    (p : ℕ) [Fact p.Prime] [CharP (ResidueField B) p]
    {π' : B} (hπ' : Irreducible π')
    (hresA : ∀ b : B, ∃ a : A, b - algebraMap A B a ∈ maximalIdeal B)
    (hadj : Algebra.adjoin A ({π'} : Set B) = ⊤)
    (hAinj : Function.Injective (algebraMap A B))
    (h0 : lowerRamificationGroup B G 0 = ⊤) (habel : ∀ x y : G, x * y = y * x)
    (k : ℕ) :
    Nat.card (G ⧸ upperRamificationGroup G π' ((k + 1 : ℕ) : ℝ))
      ∣ (Nat.card (ResidueField B) - 1) * Nat.card (ResidueField B) ^ k := by
  have h := natCard_quot_upperRamificationGroup_dvd (A := A) p hπ' hresA hadj hAinj h0 habel (k + 1)
  simpa using h

end Main

end ABC3.Found.PGC
