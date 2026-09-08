import ABC3.Found.PGC.GainedTowerDescent
import ABC3.Found.PGC.MinpolyOrbitSplit
import Mathlib.Data.Nat.Choose.Dvd

/-!
# [pGC] ★★★★降下の値段の模型を**体の層へ橋渡しする** —— `gainedLoss` は実際の損失を上から押さえる

`Found/PGC/GainedTowerDescent.lean` は `ℤ` 上の漸化式 `gainedLoss` を立て、
`(p−1)²Λ_k ≤ (p^{k+1}−p)e`(`gainedLoss_fits`)まで閉じた。そこに残っていたのは

> ★**降下の値段の模型そのものは Lean の外にある。** `gainedLoss` は `ℤ` 上の `def` であり、
> 「この漸化式が実際の降下の損失を上から押さえる」ことは測定 1〜3 の**議論**であって
> Lean の証明ではない。

の 1 点である。★★**本ファイルはこの橋を架ける。** 最終の出口は

```
exists_norm_sub_algebraMap_le_prod_axDecay :
  ∃ y ∈ F,  ‖x − y‖ ≤ (∏_{j∈[1,k+1]} axDecay p j) · ‖τ x − x‖
```

すなわち★★**体の元 `x` に対する `AxLemmaGraded` の形そのもの**である。

## ★★★測定 0 —— 配られた字面の 3 本のうち **1 本目は木に既に在った**

持ち場は「実装者が挙げた依拠(★木にも mathlib にも無いと実測済み)」として 3 本を挙げていた。
★★**1 本目は偽である。** 測定:

```
grep -n "norm_sub_digit_zero_eq_zpow_mul" lean/ABC3/Found/PGC/CyclicJumpNorm.lean
```

`CyclicJumpNorm.norm_sub_digit_zero_eq_zpow_mul`(542 行)が

```
‖digitSum n π a − algebraMap K L (a 0)‖ = ‖π‖^{−i} · ‖σ (digitSum n π a) − digitSum n π a‖
```

を**等式で**持っている。これは持ち場の 1 本目

  `max_{y∈E} v_M(x−y) = v_M((σ−1)x) − t`

をノルム言語で書いたもの**そのもの**であり、`exists_norm_sub_algebraMap_eq_zpow_mul`
(575 行)は「どの `c ∈ K` もこれより近くならない」まで付けている。
★`GainedTowerDescent` の docstring が「木にも mathlib にも無い」と書いていたのは
Serre V §3 Lemme 4 / IV §1 Prop 4 / III §6 Prop 13 の 3 本であって、
★**1 層の降下の等式はそこに入っていない**(持ち場の転記が 1 本ずれている)。

## ★★★★測定 1 —— 本当に欠けていたのは `(τ^p − 1) ≡ (τ−1)^p (mod p)` の方だった

`GainedTowerDescent` の測定 1(得の補題)は docstring の議論のままだった。★本ファイルの §1・§2 が
これを**分岐・付値・Galois の語彙が 1 語も出ない抽象核**として形式化する:

```
norm_iterate_prime_sub_self_le :
  τ : L →+ L,  ‖τ z − z‖ ≤ θ‖z‖ (∀z),  0 ≤ θ ≤ 1,  p 素数
    ⟹  ‖τ^[p] x − x‖ ≤ max(θ^{p−1}, ‖(p:L)‖) · ‖τ x − x‖
```

★仮定に要るのは「超距離ノルム体」「`τ` が**加法的**」「収縮率 `θ`」だけである。
★`τ` の乗法性も等長性も全単射性も使わない(`AddMonoidHom` で足りる)。
★★証明は `τ^[n] x = Σ_i C(n,i)·Δ^[i] x`(`iterate_eq_sum_choose_smul`、Newton の前進差分)と
`p ∣ C(p,m)`(`0<m<p`)の 2 つだけ。

★塔に沿って繰り返した形が `norm_iterate_pow_sub_self_le`:

```
‖τ^[p^k] x − x‖ ≤ (∏_{j<k} θ_j^{p−1}) · ‖τ x − x‖    (θ_j は τ^{p^j} の収縮率)
```

★★ここで `max` が `θ_j^{p−1}` に潰れる条件が `hlayer : ‖(p:L)‖ ≤ θ_j^{p−1}`、
すなわち★**sharp な層の上界 `(p−1)t_j ≤ e_M`** である
(`RamificationJumpBound.norm_natCast_le_pow_of_prod_X_sub_C` がノルム言語で出している形)。

## ★★★測定 2 —— 3 本目(Serre V §3 Lemme 4、跡の像)は**要らなかった**

持ち場は「`Tr(𝔭_L^n) = 𝔭_M^{⌊(n+d)/e⌋}` が一番重い」と見ていた。★★**本ファイルは 1 度も使わない。**
理由は `GainedTowerDescent` の測定 2 が既に書いていた:

> ★欠損 `γ = v_M(p) − (p−1)t` は `(1/p)Tr` を使うから発生していたのであって、
> **降下そのものには不要**であった。

★降下は桁展開の定数項 `a_0` を取るだけ(`CyclicJumpNorm`)で、跡も射影も通らない。
★2 本目(Serre IV §1 Prop 4、`d = (p−1)(t+1)`)も**直接には使わない**: 必要なのは
その帰結 `(p−1)t ≤ e` の**ノルム版**だけで、それは `RamificationJumpBound` が
最小多項式の微分を 2 通りに測って既に出している。★★**残っていたのは 1 本もなく、
欠けていたのは「配られていなかった 4 本目」(得の補題)だった。**

## ★★本ファイルが出したもの(すべて `sorry` 0)

| 宣言 | 内容 | 抽象核か |
|---|---|---|
| `iterate_eq_sum_choose_smul` | `τ^[n]x = Σ_i C(n,i)·Δ^[i]x` | ★加群のみ |
| `norm_iterate_delta_le` | `‖Δ^[j]z‖ ≤ θ^j‖z‖` | ★ノルムのみ |
| ★★★`norm_iterate_prime_sub_self_le` | ★**得の補題**(`max(θ^{p−1}, ‖p‖)`) | ★★分岐 0 語 |
| ★★★`norm_iterate_pow_sub_self_le` | 塔に沿った得の積 | ★★分岐 0 語 |
| ★★`norm_sub_digit_zero_le_gain` | 1 層の等式 × 得 | 具体層 |
| ★★★`norm_sub_digit_zero_le_topDefect` | ★指数がちょうど `A_k = t_k − (p−1)G_{k−1}` | 具体層 |
| ★★★`norm_sub_digit_zero_le_gainedLoss` | ★**`‖x−y‖ ≤ ‖π‖^{−Λ_k}·ε`** | 具体層 |
| ★★★★`exists_norm_sub_algebraMap_le_prod_axDecay` | ★**`AxLemmaGraded` の形** | 具体層 |

## ★★指数の勘定が合っていること(★数値で検算した)

`ℚ₃(ζ₂₇)/ℚ₃(ζ₃)`: `p = 3`, `e = e_K = 2`, `k + 1 = 2`, `e_M = p²e = 18`。
下付き跳び(`v_M`)は `t_1 = 2`(`i_G(τ) = 3`)、`t_2 = 8`(`i_G(τ³) = 9`)。

* `hstepZ`: `3·2 = 6 ≤ 8` ✓、`hlayerZ`: `2·2 = 4 ≤ 6`、`2·8 = 16 ≤ 18` ✓
* `hlayer`(ノルム版): `‖3‖ = ‖π‖^{18} ≤ ‖π‖^{4}`、`≤ ‖π‖^{16}` ✓
* `Λ_2 = max(t_2, p·t_1) = max(8, 6) = 8`(`gainedLoss_two`)
* 結論の左: `‖π‖^{−8} = 3^{8/18} = 3^{4/9}`、右: `∏_{j∈[1,2]} axDecay 3 j = 3^{(1/2)(1+1/3)} = 3^{2/3}`
* `4/9 ≤ 6/9` ✓ ★**`Λ_2 = 8` は `EquivariantProjectionDescent` の実測 `L = 8` と一致する。**

## ★★★閉じていないもの(★正確に)

1. ★**`hstep` / `hlayer` / `hbreak` / `hval` / `hchar` は仮説のままである。**
   これらは「`τ^{p^j}` が下付き分岐群 `Γ_{t_{j+1}}` に入る」「層の sharp な上界」
   「`M = F(π)` が全分岐」の**ノルム言語での言い換え**であり、
   ★本ファイルは*それらから先*を閉じた。分岐群の定義からこれらを出す部分は本ファイルの外である。
2. ★**`hlayer`(`(p−1)t_j ≤ e_M`)と `hlayerZ`(`(p−1)t_j ≤ p^j e`)を別々に取っている。**
   `j ≤ k` では後者が前者を含むが、`hlayer` は塔の帰納のために `∀ j` で要るので
   (`j > k` では `p^{j}e > e_M`)、両方を仮説にした。★逸脱ではなく、仮説が 1 つ多い。
3. ★**`gainedLoss` の第 2 枝(下の塔まで降りる `max(0,A_m−t_1) + p·Λ_{m−1}`)は使っていない。**
   本ファイルが体の層で実現したのは★**第 1 枝 `A_m` の方**であり、
   `topDefect_le_gainedLoss`(`A_m ≤ Λ_m`)で `Λ` に緩めて出口を合わせている。
   ★第 2 枝を体の層で実現するには中間体の塔 `F_0 ⊂ ⋯ ⊂ F_k` を作る必要があり、
   それは `IntermediateField` の 2 層をまたぐ道(`lean-idioms` #59/#69)に当たるので本波では避けた。
   ★★**`A_m ≤ Λ_m` なので出口(`gainedLoss_fits` の予算)には影響しない。**
4. ★塔が**巡回でない**場合(`ℚ₂(ζ₁₆)`)は `τ` が 1 本で取れないので適用外(`GainedTowerDescent` と同じ)。

## 逸脱の記録(CLAUDE.md「逸脱」)

1. 原典(Serre *Corps Locaux* IV, Ax 1970)は付値言語。本ファイルはノルム言語
   (`CyclicJumpNorm` / `AxEpsilonDecay` に合わせる)。`v_M(z) = m` は `‖z‖ = ‖π‖^m`。
2. ★`τ` を `M →+ M`(加法準同型)としてしか使わない。原典は Galois 群の元だが、
   ★得の補題には加法性しか要らないので**仮定を落とした**。
   具体層で `σ : M →ₐ[F] M` と繋ぐときだけ `hσ : ∀ z, σ z = τ^[p^k] z` を要求する。
3. `GainedTowerDescent` / `CyclicJumpNorm` / `MinpolyOrbitSplit` は**読むだけ**で書き換えない。
-/

open Finset

namespace ABC3.Found.PGC

namespace GainedBridge

/-! ## §1 抽象核(その 1) —— Newton の前進差分

★ここには分岐・付値・Galois・体すら出てこない。可換加法群と加法準同型だけである。 -/

section Diff

variable {A : Type*} [AddCommGroup A]

/-- 差分作用素 `Δ_τ = τ − id`。★`τ` の乗法性は要らない。 -/
def deltaHom (τ : A →+ A) : A →+ A := τ - AddMonoidHom.id A

@[simp] theorem deltaHom_apply (τ : A →+ A) (z : A) : deltaHom τ z = τ z - z := rfl

/-- `A →+ A` を `AddMonoid.End A` として見る。

★★これが要る理由(実測): `AddMonoid.End A` は `A →+ A` の**別名**だが、
型注釈 `(τ : AddMonoid.End A)` だけではインスタンス探索が切り替わらず

```
failed to synthesize instance of type class
  HPow (L →+ L) ℕ ?m.113
```

で落ちる。★`def` の返り値型で渡すと通る。 -/
def toEnd (τ : A →+ A) : AddMonoid.End A := τ

@[simp] theorem coe_toEnd (τ : A →+ A) : ⇑(toEnd τ) = ⇑τ := rfl

/-- ★★**Newton の前進差分**: `τ^[n] x = Σ_{i≤n} C(n,i)·Δ^[i] x`。

★`(1 + Δ)^n = Σ C(n,i)Δ^i` の作用素恒等式を、`AddMonoid.End` の環構造を通さずに
`x` の上で直接帰納したもの。Pascal(`Nat.choose_succ_succ`)しか使わない。 -/
theorem iterate_eq_sum_choose_smul (τ : A →+ A) : ∀ (n : ℕ) (x : A),
    (⇑τ)^[n] x = ∑ i ∈ Finset.range (n + 1), n.choose i • (⇑(deltaHom τ))^[i] x := by
  intro n
  induction n with
  | zero => intro x; simp
  | succ n ih =>
    intro x
    have hstep : (⇑τ)^[n + 1] x = (⇑τ)^[n] x + deltaHom τ ((⇑τ)^[n] x) := by
      rw [Function.iterate_succ_apply', deltaHom_apply]; abel
    have hD : deltaHom τ (∑ i ∈ Finset.range (n + 1), n.choose i • (⇑(deltaHom τ))^[i] x)
        = ∑ i ∈ Finset.range (n + 1), n.choose i • (⇑(deltaHom τ))^[i + 1] x := by
      rw [map_sum]
      exact Finset.sum_congr rfl fun i _ => by rw [map_nsmul, Function.iterate_succ_apply']
    rw [hstep, ih x, hD]
    rw [Finset.sum_range_succ' (fun i => (n + 1).choose i • (⇑(deltaHom τ))^[i] x) (n + 1)]
    rw [Finset.sum_range_succ' (fun i => n.choose i • (⇑(deltaHom τ))^[i] x) n]
    have hpascal : ∀ i, (n + 1).choose (i + 1) • (⇑(deltaHom τ))^[i + 1] x
        = n.choose i • (⇑(deltaHom τ))^[i + 1] x
          + n.choose (i + 1) • (⇑(deltaHom τ))^[i + 1] x := fun i => by
      rw [Nat.choose_succ_succ, add_smul]
    simp only [hpascal, Finset.sum_add_distrib]
    rw [Finset.sum_range_succ (fun i => n.choose (i + 1) • (⇑(deltaHom τ))^[i + 1] x) n]
    simp [Nat.choose_succ_self]
    abel

/-- `i = 0` の項(= `x`)を落とした形。★`τ^[n]x − x` は `Δ^[1]` 以降だけで書ける。 -/
theorem iterate_sub_self_eq_sum (τ : A →+ A) (n : ℕ) (x : A) :
    (⇑τ)^[n] x - x = ∑ i ∈ Finset.range n, n.choose (i + 1) • (⇑(deltaHom τ))^[i + 1] x := by
  rw [iterate_eq_sum_choose_smul τ n x,
    Finset.sum_range_succ' (fun i => n.choose i • (⇑(deltaHom τ))^[i] x) n]
  simp

end Diff

/-! ## §2 抽象核(その 2) —— ★★得の補題

★★**「`τ` でほぼ不変」な `x` は、`τ^p` では `θ^{p−1}` だけ**得をしている**。**
分岐・付値・Galois の語彙が 1 語も出ない。 -/

section NormCore

variable {L : Type*} [NormedField L] [IsUltrametricDist L]

/-- `p ∣ c` なら `‖(c:L)‖ ≤ ‖(p:L)‖`(超距離では `‖(m:L)‖ ≤ 1`)。

★`IsUltrametricDist.norm_natCast_le_one` は section の `variable (R)` を含むので
★**明示引数が 1 つ多い**(`lean-idioms` #297 と同じ形)。`... L m` と書く。 -/
theorem norm_natCast_le_of_dvd {p c : ℕ} (hdvd : p ∣ c) : ‖(c : L)‖ ≤ ‖(p : L)‖ := by
  obtain ⟨m, rfl⟩ := hdvd
  rw [Nat.cast_mul, norm_mul]
  calc ‖(p : L)‖ * ‖(m : L)‖ ≤ ‖(p : L)‖ * 1 :=
        mul_le_mul_of_nonneg_left (IsUltrametricDist.norm_natCast_le_one L m) (norm_nonneg _)
    _ = ‖(p : L)‖ := mul_one _

omit [IsUltrametricDist L] in
/-- 収縮率 `θ` は反復で `θ^j` になる。 -/
theorem norm_iterate_delta_le {τ : L →+ L} {θ : ℝ} (hθ0 : 0 ≤ θ)
    (hD : ∀ z : L, ‖τ z - z‖ ≤ θ * ‖z‖) :
    ∀ (j : ℕ) (z : L), ‖(⇑(deltaHom τ))^[j] z‖ ≤ θ ^ j * ‖z‖ := by
  intro j
  induction j with
  | zero => intro z; simp
  | succ j ih =>
    intro z
    rw [Function.iterate_succ_apply]
    calc ‖(⇑(deltaHom τ))^[j] (deltaHom τ z)‖ ≤ θ ^ j * ‖deltaHom τ z‖ := ih _
      _ ≤ θ ^ j * (θ * ‖z‖) :=
          mul_le_mul_of_nonneg_left (by simpa using hD z) (pow_nonneg hθ0 j)
      _ = θ ^ (j + 1) * ‖z‖ := by ring

/-- ★★★★**得の補題**(`GainedTowerDescent` 測定 1 の形式化)。

`τ` が収縮率 `θ ≤ 1` の加法作用素なら

  `‖τ^[p] x − x‖ ≤ max(θ^{p−1}, ‖(p:L)‖) · ‖τ x − x‖`.

★中身は `τ^p − 1 = Δ^p + Σ_{m=1}^{p−1} C(p,m)Δ^{p−m}` と `p ∣ C(p,m)` だけ:
第 1 項が `θ^{p−1}‖Δx‖`、残りが `‖p‖·‖Δx‖` 以下で、超距離が和を `max` に潰す。

★★**`p ≥ 5` の窓を消したのはこの補題である**(旧上界はこの得を 1 度も使っていなかった)。 -/
theorem norm_iterate_prime_sub_self_le {p : ℕ} (hp : p.Prime) {τ : L →+ L} {θ : ℝ}
    (hθ0 : 0 ≤ θ) (hθ1 : θ ≤ 1) (hD : ∀ z : L, ‖τ z - z‖ ≤ θ * ‖z‖) (x : L) :
    ‖(⇑τ)^[p] x - x‖ ≤ max (θ ^ (p - 1)) ‖(p : L)‖ * ‖τ x - x‖ := by
  have hmax0 : (0 : ℝ) ≤ max (θ ^ (p - 1)) ‖(p : L)‖ :=
    le_trans (norm_nonneg _) (le_max_right _ _)
  have hC0 : (0 : ℝ) ≤ max (θ ^ (p - 1)) ‖(p : L)‖ * ‖τ x - x‖ :=
    mul_nonneg hmax0 (norm_nonneg _)
  rw [iterate_sub_self_eq_sum]
  refine IsUltrametricDist.norm_sum_le_of_forall_le_of_nonneg hC0 fun i hi => ?_
  rw [Finset.mem_range] at hi
  have hiter : (⇑(deltaHom τ))^[i + 1] x = (⇑(deltaHom τ))^[i] (τ x - x) := by
    rw [Function.iterate_succ_apply, deltaHom_apply]
  have hnorm : ‖(⇑(deltaHom τ))^[i + 1] x‖ ≤ θ ^ i * ‖τ x - x‖ := by
    rw [hiter]; exact norm_iterate_delta_le hθ0 hD i _
  rw [nsmul_eq_mul, norm_mul]
  rcases eq_or_lt_of_le (Nat.succ_le_of_lt hi) with heq | hlt
  · -- `i + 1 = p` の項。係数は `C(p,p) = 1`、`Δ^p` が `θ^{p−1}` を稼ぐ。
    have heq' : i + 1 = p := heq
    have hc : p.choose (i + 1) = 1 := by rw [heq', Nat.choose_self]
    have hip : θ ^ i = θ ^ (p - 1) := by rw [show i = p - 1 by omega]
    rw [hc, Nat.cast_one, norm_one, one_mul]
    calc ‖(⇑(deltaHom τ))^[i + 1] x‖ ≤ θ ^ i * ‖τ x - x‖ := hnorm
      _ = θ ^ (p - 1) * ‖τ x - x‖ := by rw [hip]
      _ ≤ max (θ ^ (p - 1)) ‖(p : L)‖ * ‖τ x - x‖ :=
          mul_le_mul_of_nonneg_right (le_max_left _ _) (norm_nonneg _)
  · -- `1 ≤ i + 1 < p` の項。★係数が `p` で割れる。
    have hdvd : p ∣ p.choose (i + 1) := hp.dvd_choose_self (Nat.succ_ne_zero i) hlt
    have h1 : ‖((p.choose (i + 1) : ℕ) : L)‖ ≤ ‖(p : L)‖ := norm_natCast_le_of_dvd hdvd
    have h2 : ‖(⇑(deltaHom τ))^[i + 1] x‖ ≤ ‖τ x - x‖ := by
      refine le_trans hnorm ?_
      have hle : θ ^ i ≤ 1 := pow_le_one₀ hθ0 hθ1
      nlinarith [norm_nonneg (τ x - x), pow_nonneg hθ0 i]
    calc ‖((p.choose (i + 1) : ℕ) : L)‖ * ‖(⇑(deltaHom τ))^[i + 1] x‖
        ≤ ‖(p : L)‖ * ‖τ x - x‖ := mul_le_mul h1 h2 (norm_nonneg _) (norm_nonneg _)
      _ ≤ max (θ ^ (p - 1)) ‖(p : L)‖ * ‖τ x - x‖ :=
          mul_le_mul_of_nonneg_right (le_max_right _ _) (norm_nonneg _)

/-- ★★★塔に沿って得の補題を繰り返した形。

`θ j` は `τ^{p^j}` の収縮率(具体層では `‖π‖^{t_{j+1}}`)。`hlayer` は sharp な層の上界
`‖(p:L)‖ ≤ θ_j^{p−1}` で、これが `max` を `θ_j^{p−1}` に潰す。

  `‖τ^[p^k] x − x‖ ≤ (∏_{j<k} θ_j^{p−1}) · ‖τ x − x‖`

★分岐・付値・Galois の語彙が 1 語も出ない。 -/
theorem norm_iterate_pow_sub_self_le {p : ℕ} (hp : p.Prime) {τ : L →+ L} {θ : ℕ → ℝ}
    (hθ0 : ∀ j, 0 ≤ θ j) (hθ1 : ∀ j, θ j ≤ 1)
    (hstep : ∀ (j : ℕ) (z : L), ‖(⇑τ)^[p ^ j] z - z‖ ≤ θ j * ‖z‖)
    (hlayer : ∀ j, ‖(p : L)‖ ≤ θ j ^ (p - 1)) :
    ∀ (k : ℕ) (x : L),
      ‖(⇑τ)^[p ^ k] x - x‖ ≤ (∏ j ∈ Finset.range k, θ j ^ (p - 1)) * ‖τ x - x‖ := by
  intro k
  induction k with
  | zero => intro x; simp
  | succ k ih =>
    intro x
    have hg : ⇑((toEnd τ) ^ (p ^ k)) = (⇑τ)^[p ^ k] :=
      AddMonoid.End.coe_pow L (toEnd τ) (p ^ k)
    have hiter : (⇑τ)^[p ^ (k + 1)] x = (⇑((toEnd τ) ^ (p ^ k)))^[p] x := by
      rw [hg, ← Function.iterate_mul, ← pow_succ]
    have hD : ∀ z : L, ‖((toEnd τ) ^ (p ^ k)) z - z‖ ≤ θ k * ‖z‖ := by
      intro z; rw [hg]; exact hstep k z
    have h1 := norm_iterate_prime_sub_self_le hp (τ := (toEnd τ) ^ (p ^ k))
      (hθ0 k) (hθ1 k) hD x
    rw [max_eq_left (hlayer k)] at h1
    rw [hiter]
    have h2 : ‖((toEnd τ) ^ (p ^ k)) x - x‖
        ≤ (∏ j ∈ Finset.range k, θ j ^ (p - 1)) * ‖τ x - x‖ := by
      rw [hg]; exact ih x
    have hpow0 : (0 : ℝ) ≤ θ k ^ (p - 1) := pow_nonneg (hθ0 k) _
    calc ‖(⇑((toEnd τ) ^ (p ^ k)))^[p] x - x‖
        ≤ θ k ^ (p - 1) * ‖((toEnd τ) ^ (p ^ k)) x - x‖ := h1
      _ ≤ θ k ^ (p - 1) * ((∏ j ∈ Finset.range k, θ j ^ (p - 1)) * ‖τ x - x‖) :=
          mul_le_mul_of_nonneg_left h2 hpow0
      _ = (∏ j ∈ Finset.range (k + 1), θ j ^ (p - 1)) * ‖τ x - x‖ := by
          rw [Finset.prod_range_succ]; ring

end NormCore

/-! ## §3 体の層への橋 —— 1 層の降下の**等式** × 塔に沿った**得**

★1 層の降下の等式は木に既に在る(`CyclicJumpNorm.norm_sub_digit_zero_eq_zpow_mul`)。
本節はそこに §2 の得を掛けるだけである。 -/

section Bridge

variable {F M : Type*} [Field F] [NormedField M] [IsUltrametricDist M] [Algebra F M]

/-- ★★1 層の降下(等式)に塔の得を掛けた形。

`σ = τ^{p^k}` は層 `M/F` の生成元、`τ` は塔 `M/K` の生成元。跳び `i` は `v_M` の目盛りで
`‖σπ − π‖ = ‖π‖^{i+1}`。結論は

  `‖x − a_0‖ ≤ ‖π‖^{−i} · (∏_{j<k} θ_j^{p−1}) · ‖τ x − x‖`

★左辺は `d(x, F)` そのもの(`CyclicJumpNorm` が最良近似であることまで示している)。 -/
theorem norm_sub_digit_zero_le_gain {p i k : ℕ} (hp : p.Prime) {π : M} {a : ℕ → F}
    (τ : M →+ M) (σ : M →ₐ[F] M) {θ : ℕ → ℝ}
    (hσ : ∀ z : M, σ z = (⇑τ)^[p ^ k] z)
    (hθ0 : ∀ j, 0 ≤ θ j) (hθ1 : ∀ j, θ j ≤ 1)
    (hstep : ∀ (j : ℕ) (z : M), ‖(⇑τ)^[p ^ j] z - z‖ ≤ θ j * ‖z‖)
    (hlayer : ∀ j, ‖(p : M)‖ ≤ θ j ^ (p - 1))
    (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1) (hi : 0 < i)
    (hbreak : ‖σ π - π‖ = ‖π‖ ^ (i + 1))
    (hchar : ∀ j, 0 < j → j < p → ‖(j : M)‖ = 1)
    (hval : ∀ c : F, c ≠ 0 → ∃ m : ℤ, ‖algebraMap F M c‖ = ‖π‖ ^ ((p : ℤ) * m)) :
    ‖digitSum p π a - algebraMap F M (a 0)‖
      ≤ ‖π‖ ^ (-(i : ℤ)) * (∏ j ∈ Finset.range k, θ j ^ (p - 1))
          * ‖τ (digitSum p π a) - digitSum p π a‖ := by
  have heq := norm_sub_digit_zero_eq_zpow_mul (a := a) σ hp.pos hπ0 hπ1 hi hbreak hchar hval
  have hchain := norm_iterate_pow_sub_self_le hp hθ0 hθ1 hstep hlayer k (digitSum p π a)
  rw [heq, hσ (digitSum p π a)]
  have hz0 : (0 : ℝ) ≤ ‖π‖ ^ (-(i : ℤ)) := le_of_lt (zpow_pos hπ0 _)
  calc ‖π‖ ^ (-(i : ℤ)) * ‖(⇑τ)^[p ^ k] (digitSum p π a) - digitSum p π a‖
      ≤ ‖π‖ ^ (-(i : ℤ)) * ((∏ j ∈ Finset.range k, θ j ^ (p - 1))
          * ‖τ (digitSum p π a) - digitSum p π a‖) := mul_le_mul_of_nonneg_left hchain hz0
    _ = _ := by ring

end Bridge

/-! ## §4 指数の勘定 —— `‖π‖^{−i}·∏θ_j^{p−1} = ‖π‖^{−A_k}`

`A_k = t_k − (p−1)·G_{k−1}` は `gainedLoss` の第 1 枝そのものである。 -/

section Exponent

/-- 抽象核(実数のみ): `∏_{j<k} b^{c_j} = b^{Σ c_j}`。 -/
theorem prod_zpow_range (b : ℝ) (hb : 0 < b) (c : ℕ → ℤ) : ∀ k : ℕ,
    (∏ j ∈ Finset.range k, b ^ (c j)) = b ^ (∑ j ∈ Finset.range k, c j) := by
  intro k
  induction k with
  | zero => simp
  | succ k ih => rw [Finset.prod_range_succ, Finset.sum_range_succ, ih, zpow_add₀ (ne_of_gt hb)]

/-- `jumpSum` は `Finset.range` 上の和である。 -/
theorem jumpSum_eq_sum_range (t : ℕ → ℤ) : ∀ k : ℕ,
    GainedDescent.jumpSum t k = ∑ j ∈ Finset.range k, t (j + 1) := by
  intro k
  induction k with
  | zero => simp [GainedDescent.jumpSum]
  | succ k ih => rw [GainedDescent.jumpSum, ih, Finset.sum_range_succ]

/-- ★`θ_j = b^{t_{j+1}}` のとき、得の積は `b^{(p−1)·G_k}` である。 -/
theorem prod_theta_eq {b : ℝ} (hb : 0 < b) {p : ℕ} (hp : 1 ≤ p) (t : ℕ → ℤ) (k : ℕ) :
    (∏ j ∈ Finset.range k, (b ^ (t (j + 1))) ^ (p - 1))
      = b ^ (((p : ℤ) - 1) * GainedDescent.jumpSum t k) := by
  have hcast : ((p - 1 : ℕ) : ℤ) = (p : ℤ) - 1 := by
    push_cast [Nat.cast_sub hp]; ring
  have hterm : ∀ j : ℕ, (b ^ (t (j + 1))) ^ (p - 1) = b ^ (((p : ℤ) - 1) * t (j + 1)) := by
    intro j
    rw [← zpow_natCast (b ^ (t (j + 1))) (p - 1), ← zpow_mul, hcast]
    congr 1
    ring
  simp only [hterm]
  rw [prod_zpow_range b hb (fun j => ((p : ℤ) - 1) * t (j + 1)) k, jumpSum_eq_sum_range,
    Finset.mul_sum]

/-- ★`gainedLoss` の第 1 枝: `A_{m+1} ≤ Λ_{m+1}`。 -/
theorem topDefect_le_gainedLoss (p : ℤ) (t : ℕ → ℤ) (m : ℕ) :
    t (m + 1) - (p - 1) * GainedDescent.jumpSum t m ≤ GainedDescent.gainedLoss p t (m + 1) := by
  rw [GainedDescent.gainedLoss]
  exact le_max_left _ _

/-- `p·t_j ≤ t_{j+1}` のもとで `Λ_{m+1} ≥ 0`。 -/
theorem gainedLoss_nonneg {p : ℤ} {t : ℕ → ℤ} (ht1 : 0 ≤ t 1)
    (hstep : ∀ j, p * t (j + 1) ≤ t (j + 2)) (m : ℕ) :
    0 ≤ GainedDescent.gainedLoss p t (m + 1) := by
  have h1 := GainedDescent.topDefect_nonneg ht1 hstep m
  have h2 := topDefect_le_gainedLoss p t m
  omega

/-- 抽象核(実数のみ): `b = p^{−1/E}` なら `b^{−m} = p^{m/E}`。 -/
theorem zpow_eq_rpow_div {b : ℝ} {p : ℕ} {E : ℝ} (hb : b = (p : ℝ) ^ (-(1 : ℝ) / E)) (m : ℤ) :
    b ^ (-m) = (p : ℝ) ^ ((m : ℝ) / E) := by
  rw [hb, ← Real.rpow_intCast ((p : ℝ) ^ (-(1 : ℝ) / E)) (-m),
    ← Real.rpow_mul (Nat.cast_nonneg p)]
  congr 1
  push_cast
  ring

end Exponent

/-! ## §5 組み立て —— ★★★体の層で `gainedLoss` が損失を押さえる -/

section Assemble

variable {F M : Type*} [Field F] [NormedField M] [IsUltrametricDist M] [Algebra F M]

/-- ★★★**指数がちょうど `A_k = t_k − (p−1)·G_{k−1}` になる**。

`θ_j = ‖π‖^{t_{j+1}}` と置いた §3。★これが `GainedTowerDescent` の測定 3 の段 1・2
(「上の層で止まる」枝)を体の層で実現した形である。 -/
theorem norm_sub_digit_zero_le_topDefect {p i k : ℕ} (hp : p.Prime) {π : M} {a : ℕ → F}
    (τ : M →+ M) (σ : M →ₐ[F] M) (t : ℕ → ℤ)
    (hσ : ∀ z : M, σ z = (⇑τ)^[p ^ k] z) (ht0 : ∀ j, 0 ≤ t (j + 1))
    (hstep : ∀ (j : ℕ) (z : M), ‖(⇑τ)^[p ^ j] z - z‖ ≤ ‖π‖ ^ (t (j + 1)) * ‖z‖)
    (hlayer : ∀ j, ‖(p : M)‖ ≤ (‖π‖ ^ (t (j + 1))) ^ (p - 1))
    (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1) (hi : 0 < i) (hti : (i : ℤ) = t (k + 1))
    (hbreak : ‖σ π - π‖ = ‖π‖ ^ (i + 1))
    (hchar : ∀ j, 0 < j → j < p → ‖(j : M)‖ = 1)
    (hval : ∀ c : F, c ≠ 0 → ∃ m : ℤ, ‖algebraMap F M c‖ = ‖π‖ ^ ((p : ℤ) * m)) :
    ‖digitSum p π a - algebraMap F M (a 0)‖
      ≤ ‖π‖ ^ (-(t (k + 1) - ((p : ℤ) - 1) * GainedDescent.jumpSum t k))
          * ‖τ (digitSum p π a) - digitSum p π a‖ := by
  have h := norm_sub_digit_zero_le_gain (a := a) (θ := fun j => ‖π‖ ^ (t (j + 1)))
    hp τ σ hσ (fun j => le_of_lt (zpow_pos hπ0 _))
    (fun j => zpow_le_one₀ hπ0 (le_of_lt hπ1) (ht0 j)) hstep hlayer hπ0 hπ1 hi hbreak hchar hval
  refine le_trans h (le_of_eq ?_)
  rw [prod_theta_eq hπ0 hp.one_lt.le t k, ← zpow_add₀ (ne_of_gt hπ0), hti]
  congr 1
  ring

/-- ★★★★**橋の本体**: 体の層の損失が `ℤ` の模型 `gainedLoss` で上から押さえられる。

  `d(x, F) ≤ ‖π‖^{−Λ_{k+1}} · ‖τ x − x‖`,   `Λ = gainedLoss p t (k+1)`.

★`A ≤ Λ` は `topDefect_le_gainedLoss`、`‖π‖ < 1` なので指数を大きくすると値は小さくなる。 -/
theorem norm_sub_digit_zero_le_gainedLoss {p i k : ℕ} (hp : p.Prime) {π : M} {a : ℕ → F}
    (τ : M →+ M) (σ : M →ₐ[F] M) (t : ℕ → ℤ)
    (hσ : ∀ z : M, σ z = (⇑τ)^[p ^ k] z) (ht0 : ∀ j, 0 ≤ t (j + 1))
    (hstep : ∀ (j : ℕ) (z : M), ‖(⇑τ)^[p ^ j] z - z‖ ≤ ‖π‖ ^ (t (j + 1)) * ‖z‖)
    (hlayer : ∀ j, ‖(p : M)‖ ≤ (‖π‖ ^ (t (j + 1))) ^ (p - 1))
    (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1) (hi : 0 < i) (hti : (i : ℤ) = t (k + 1))
    (hbreak : ‖σ π - π‖ = ‖π‖ ^ (i + 1))
    (hchar : ∀ j, 0 < j → j < p → ‖(j : M)‖ = 1)
    (hval : ∀ c : F, c ≠ 0 → ∃ m : ℤ, ‖algebraMap F M c‖ = ‖π‖ ^ ((p : ℤ) * m)) :
    ‖digitSum p π a - algebraMap F M (a 0)‖
      ≤ ‖π‖ ^ (-(GainedDescent.gainedLoss (p : ℤ) t (k + 1)))
          * ‖τ (digitSum p π a) - digitSum p π a‖ := by
  refine le_trans (norm_sub_digit_zero_le_topDefect hp τ σ t hσ ht0 hstep hlayer hπ0 hπ1 hi hti
    hbreak hchar hval) (mul_le_mul_of_nonneg_right ?_ (norm_nonneg _))
  refine zpow_le_zpow_right_of_le_one₀ hπ0 (le_of_lt hπ1) ?_
  have := topDefect_le_gainedLoss (p : ℤ) t k
  omega

/-- ★★桁展開を仮定しない形。`M = F(π)`(`Algebra.adjoin F {π} = ⊤`, `[M:F] = p`)なら
`M` のどの元 `x` にも `y ∈ F` が取れる。 -/
theorem exists_norm_sub_algebraMap_le_gainedLoss {p i k : ℕ} (hp : p.Prime) {π : M}
    (τ : M →+ M) (σ : M →ₐ[F] M) (t : ℕ → ℤ)
    (hσ : ∀ z : M, σ z = (⇑τ)^[p ^ k] z) (ht0 : ∀ j, 0 ≤ t (j + 1))
    (hstep : ∀ (j : ℕ) (z : M), ‖(⇑τ)^[p ^ j] z - z‖ ≤ ‖π‖ ^ (t (j + 1)) * ‖z‖)
    (hlayer : ∀ j, ‖(p : M)‖ ≤ (‖π‖ ^ (t (j + 1))) ^ (p - 1))
    (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1) (hi : 0 < i) (hti : (i : ℤ) = t (k + 1))
    (hbreak : ‖σ π - π‖ = ‖π‖ ^ (i + 1))
    (hchar : ∀ j, 0 < j → j < p → ‖(j : M)‖ = 1)
    (hval : ∀ c : F, c ≠ 0 → ∃ m : ℤ, ‖algebraMap F M c‖ = ‖π‖ ^ ((p : ℤ) * m))
    (hdeg : (minpoly F π).natDegree = p)
    (htop : Algebra.adjoin F ({π} : Set M) = ⊤) (x : M) :
    ∃ y : F, ‖x - algebraMap F M y‖
      ≤ ‖π‖ ^ (-(GainedDescent.gainedLoss (p : ℤ) t (k + 1))) * ‖τ x - x‖ := by
  have hint : IsIntegral F π := by
    by_contra hc
    rw [minpoly.eq_zero hc, Polynomial.natDegree_zero] at hdeg
    exact hp.ne_zero hdeg.symm
  obtain ⟨a, rfl⟩ := exists_digitSum_of_adjoin_eq_top hint hdeg htop x
  exact ⟨a 0, norm_sub_digit_zero_le_gainedLoss hp τ σ t hσ ht0 hstep hlayer hπ0 hπ1 hi hti
    hbreak hchar hval⟩

/-- ★★★★**出口 —— `AxLemmaGraded` の形**。

`e_M = p^{k+1}·e`(`‖π‖ = p^{−1/e_M}`)のとき

  `∃ y ∈ F,  ‖x − y‖ ≤ (∏_{j∈[1,k+1]} axDecay p j) · ‖τ x − x‖`.

★★これは `GainedTowerDescent.rpow_le_prod_axDecay`(`ℤ` の側)と
本ファイルの橋(体の側)を繋いだもので、★**降下の値段の模型が実際の損失を上から押さえ、
しかもその値が `axDecay` の積の予算に入る**ことを 1 本で述べている。

★数値検算(`ℚ₃(ζ₂₇)/ℚ₃(ζ₃)`, `p=3, e=2, k=1, t=(2,8)`): 左 `3^{4/9}` ≤ 右 `3^{2/3}`。 -/
theorem exists_norm_sub_algebraMap_le_prod_axDecay {p : ℕ} [Fact p.Prime] {e i k : ℕ} {π : M}
    (τ : M →+ M) (σ : M →ₐ[F] M) (t : ℕ → ℤ) (he : 0 < e)
    (hσ : ∀ z : M, σ z = (⇑τ)^[p ^ k] z) (ht0 : ∀ j, 0 ≤ t (j + 1))
    (hstep : ∀ (j : ℕ) (z : M), ‖(⇑τ)^[p ^ j] z - z‖ ≤ ‖π‖ ^ (t (j + 1)) * ‖z‖)
    (hlayer : ∀ j, ‖(p : M)‖ ≤ (‖π‖ ^ (t (j + 1))) ^ (p - 1))
    (hstepZ : ∀ j, (p : ℤ) * t (j + 1) ≤ t (j + 2))
    (hlayerZ : ∀ j, ((p : ℤ) - 1) * t (j + 1) ≤ (p : ℤ) ^ (j + 1) * (e : ℤ))
    (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1) (hi : 0 < i) (hti : (i : ℤ) = t (k + 1))
    (hπE : ‖π‖ = (p : ℝ) ^ (-(1 : ℝ) / ((p ^ (k + 1) * e : ℕ) : ℝ)))
    (hbreak : ‖σ π - π‖ = ‖π‖ ^ (i + 1))
    (hchar : ∀ j, 0 < j → j < p → ‖(j : M)‖ = 1)
    (hval : ∀ c : F, c ≠ 0 → ∃ m : ℤ, ‖algebraMap F M c‖ = ‖π‖ ^ ((p : ℤ) * m))
    (hdeg : (minpoly F π).natDegree = p)
    (htop : Algebra.adjoin F ({π} : Set M) = ⊤) (x : M) :
    ∃ y : F, ‖x - algebraMap F M y‖
      ≤ (∏ j ∈ Finset.Icc 1 (k + 1), axDecay p j) * ‖τ x - x‖ := by
  obtain ⟨y, hy⟩ := exists_norm_sub_algebraMap_le_gainedLoss (Fact.out : p.Prime) τ σ t hσ ht0
    hstep hlayer hπ0 hπ1 hi hti hbreak hchar hval hdeg htop x
  refine ⟨y, le_trans hy (mul_le_mul_of_nonneg_right ?_ (norm_nonneg _))⟩
  have hΛ0 : 0 ≤ GainedDescent.gainedLoss (p : ℤ) t (k + 1) :=
    gainedLoss_nonneg (by simpa using ht0 0) hstepZ k
  have hJ : ((GainedDescent.gainedLoss (p : ℤ) t (k + 1)).toNat : ℤ)
      = GainedDescent.gainedLoss (p : ℤ) t (k + 1) := Int.toNat_of_nonneg hΛ0
  have hmain := GainedDescent.rpow_le_prod_axDecay (p := p) (e := e) (k := k + 1)
    (J := (GainedDescent.gainedLoss (p : ℤ) t (k + 1)).toNat) (t := t) he ht0 hstepZ hlayerZ
    (le_of_eq hJ)
  have hcast : ((GainedDescent.gainedLoss (p : ℤ) t (k + 1) : ℤ) : ℝ)
      = (((GainedDescent.gainedLoss (p : ℤ) t (k + 1)).toNat : ℕ) : ℝ) := by
    exact_mod_cast hJ.symm
  rw [zpow_eq_rpow_div hπE, hcast]
  exact hmain

end Assemble

/-! ## §6 `.src`(原典の対応箇所)

★得の補題は原典が畳んだ箇所(`(τ^p−1) ≡ (τ−1)^p mod p` は「formally」で畳まれる型の段)。
`CyclicJumpNorm` / `RamificationJumpBound` と同じ項目を指す。 -/

def norm_iterate_prime_sub_self_le.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def norm_iterate_pow_sub_self_le.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def norm_sub_digit_zero_le_gain.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def norm_sub_digit_zero_le_gainedLoss.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def exists_norm_sub_algebraMap_le_prod_axDecay.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-! ## §7 使っている公理の一覧 -/

#print axioms iterate_eq_sum_choose_smul
#print axioms iterate_sub_self_eq_sum
#print axioms norm_natCast_le_of_dvd
#print axioms norm_iterate_delta_le
#print axioms norm_iterate_prime_sub_self_le
#print axioms norm_iterate_pow_sub_self_le
#print axioms norm_sub_digit_zero_le_gain
#print axioms prod_zpow_range
#print axioms jumpSum_eq_sum_range
#print axioms prod_theta_eq
#print axioms topDefect_le_gainedLoss
#print axioms gainedLoss_nonneg
#print axioms zpow_eq_rpow_div
#print axioms norm_sub_digit_zero_le_topDefect
#print axioms norm_sub_digit_zero_le_gainedLoss
#print axioms exists_norm_sub_algebraMap_le_gainedLoss
#print axioms exists_norm_sub_algebraMap_le_prod_axDecay

end GainedBridge

end ABC3.Found.PGC
