import ABC3.Found.PGC.ResidueCardLowerBound
import ABC3.Found.PGC.InertiaKummerBound

/-!
# `q` は `Γ_K` から群論的に決まる(経路 C のノード G)

[pGC] Proposition 1.2 への**経路 C**(`ResearchPaper/pgc-goal.md`)の合流点。
絶対 Galois 群 `Γ_K` の**位相群としての同型類だけ**から剰余体の位数
`q = #𝓀[K]` を読み取る。

## 何を合成するか

不変量は `N_n(Γ) := #Hom_cont(Γ, ℤ/n)`
(`Found/PGC/ContinuousHomCount.lean::contHomCard`)である。経路 C の主張 (C-q)
「`p ∤ n` のとき `N_n(Γ_K) = n · gcd(n, q−1)`」は、両側から挟んで得られている:

* **下界 (G1)** `Found/PGC/ResidueCardLowerBound.lean::contHomCard_absGal_of_dvd` ——
  `n ∣ q−1` なら `N_n(Γ_K) = n²`(第 1030)。
* **上界 (G2)** `Found/PGC/InertiaKummerBound.lean::contHomCard_absGal_le` ——
  `p ∤ n` なら `N_n(Γ_K) ≤ n · gcd(n, q−1)`(第 1037)。
* **不変性 (A)** `Found/PGC/ContinuousHomCount.lean::contHomCard_congr` ——
  `N_n` は位相群の同型で不変。

この 3 つから、集合

  `S(Γ) := { n : n ≠ 0, p ∤ n, N_n(Γ) = n² }`   (`residueCardSpec`)

が `q − 1` を**最大元として持つ**(`isGreatest_residueCardSpec`)。`S` は `Γ` と `p`
だけで書けているので、`Γ_K ≃ₜ* Γ_{K'}` なら `S(Γ_K) = S(Γ_{K'})`、最大元は一意だから
`q − 1 = q' − 1`、`q, q' ≥ 2` より `q = q'`。

★下界と上界の役割分担に注意。上界だけでは `S` の元が `q−1` を超えないことしか
言えず、`q−1` 自身が `S` に属することは下界が担う。両方が要る。

## 上界の側の細部

`n ∈ S` から `n ≤ q−1` を出すのに、途中で `n ∣ q−1` まで言える
(`dvd_residueCard_sub_one_of_mem_spec`)。(G2) と `N_n = n²` から
`n · n ≤ n · gcd(n, q−1)`、`n ≠ 0` で約して `n ≤ gcd(n, q−1)`。逆向きは
`gcd(n, q−1) ∣ n` から常に成り立つので `gcd(n, q−1) = n`、すなわち `n ∣ q−1`。
つまり `S(Γ_K)` はちょうど `q−1` の約数のうち `0` でないものの全体である
(下界と合わせれば等号)。

## 惰性群を名指ししていないこと(設計上の記録)

本ファイルは `inertia` にも `Found/PGC/InertiaIdentification.lean` の
`absGalQuotKerEquivUnramifiedGal` にも触れていない。それらは `Interface` の
`SubgroupCorrespondence` / `ResidueCardinality` を経由するため、経路 C の出口である
Corollary 1.3 と循環する。経路 C は `contHomCard` という**位相群だけで書ける**
不変量を通すことで、その循環を避けている。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC ABC3.Found
open scoped NormedField Valued

variable {p : ℕ} [Fact p.Prime]

/-! ## 群論的な特徴づけの集合 -/

/-- **`q − 1` を切り出す群論的な集合** `S(Γ_K) = { n ≠ 0 : p ∤ n, N_n(Γ_K) = n² }`。

`Γ_K`(位相群)と `p` だけで書けていることが要点である。体 `K` は現れない。 -/
def residueCardSpec (K : PAdicLocalField p) : Set ℕ :=
  {n : ℕ | n ≠ 0 ∧ ¬ p ∣ n ∧ contHomCard K.absGal n = n * n}

/-! ## `q − 1` についての基本 -/

/-- `q ≥ 2` なので `q − 1 ≠ 0`。 -/
theorem residueCard_sub_one_ne_zero (K : PAdicLocalField p) :
    Nat.card 𝓀[K.carrier] - 1 ≠ 0 := by
  have h1 : 1 < Fintype.card 𝓀[K.carrier] := Fintype.one_lt_card
  rw [Nat.card_eq_fintype_card]
  omega

/-- `p ∤ q − 1`(`q = p^f`, `f ≥ 1`)の `Nat.card` 版。中身は
`Found/PGC/PrimeToPTorsion.lean::not_dvd_residueCard_sub_one`。 -/
theorem not_dvd_residueCard_sub_one' (K : PAdicLocalField p) :
    ¬ p ∣ (Nat.card 𝓀[K.carrier] - 1) := by
  rw [Nat.card_eq_fintype_card]
  exact not_dvd_residueCard_sub_one K

/-! ## `q − 1 ∈ S(Γ_K)` —— 下界 (G1) の側 -/

/-- `q − 1` は `S(Γ_K)` に属する。`n := q−1` に (G1) を `dvd_rfl` で当てるだけ。 -/
theorem residueCard_sub_one_mem_spec (K : PAdicLocalField p) :
    Nat.card 𝓀[K.carrier] - 1 ∈ residueCardSpec K :=
  ⟨residueCard_sub_one_ne_zero K, not_dvd_residueCard_sub_one' K,
    contHomCard_absGal_of_dvd K (not_dvd_residueCard_sub_one' K)
      (residueCard_sub_one_ne_zero K) dvd_rfl⟩

/-! ## `S(Γ_K)` の元は `q − 1` を割る —— 上界 (G2) の側 -/

/-- `n ∈ S(Γ_K)` なら `n ∣ q − 1`。

(G2) `N_n ≤ n · gcd(n, q−1)` と `N_n = n²` から `n ≤ gcd(n, q−1)`、逆向きは常に
成り立つので `gcd(n, q−1) = n`。 -/
theorem dvd_residueCard_sub_one_of_mem_spec (K : PAdicLocalField p) {n : ℕ}
    (hn : n ∈ residueCardSpec K) : n ∣ Nat.card 𝓀[K.carrier] - 1 := by
  obtain ⟨hn0, hnp, hcard⟩ := hn
  have h := contHomCard_absGal_le K hnp
  rw [hcard] at h
  have h2 : n ≤ Nat.gcd n (Nat.card 𝓀[K.carrier] - 1) :=
    Nat.le_of_mul_le_mul_left h (Nat.pos_of_ne_zero hn0)
  have h3 : Nat.gcd n (Nat.card 𝓀[K.carrier] - 1) = n :=
    le_antisymm (Nat.le_of_dvd (Nat.pos_of_ne_zero hn0) (Nat.gcd_dvd_left _ _)) h2
  exact h3 ▸ Nat.gcd_dvd_right _ _

/-- `n ∈ S(Γ_K)` なら `n ≤ q − 1`。 -/
theorem le_residueCard_sub_one_of_mem_spec (K : PAdicLocalField p) {n : ℕ}
    (hn : n ∈ residueCardSpec K) : n ≤ Nat.card 𝓀[K.carrier] - 1 :=
  Nat.le_of_dvd (Nat.pos_of_ne_zero (residueCard_sub_one_ne_zero K))
    (dvd_residueCard_sub_one_of_mem_spec K hn)

/-- **★★★★★★★★★★★★★★★★★★(G) `q − 1 = max { n : p ∤ n, N_n(Γ_K) = n² }`**。

経路 C の (C-q) の帰結。下界 (G1) が「属すること」を、上界 (G2) が「最大であること」を
担う。左辺は体 `K` の不変量、右辺は位相群 `Γ_K` と `p` だけで書ける量である。 -/
theorem isGreatest_residueCardSpec (K : PAdicLocalField p) :
    IsGreatest (residueCardSpec K) (Nat.card 𝓀[K.carrier] - 1) :=
  ⟨residueCard_sub_one_mem_spec K, fun _ hn => le_residueCard_sub_one_of_mem_spec K hn⟩

/-- **`S(Γ_K)` の正体** —— ちょうど「`q − 1` の `0` でない約数」の全体。

`⊆` が上界 (G2)、`⊇` が下界 (G1) である。`p ∤ n` は `n ∣ q−1` から従う
(`Found/PGC/ResidueCardLowerBound.lean::not_dvd_of_dvd_residueCard_sub_one`)ので、
右辺に `p` は現れない。 -/
theorem residueCardSpec_eq (K : PAdicLocalField p) :
    residueCardSpec K = {n : ℕ | n ≠ 0 ∧ n ∣ Nat.card 𝓀[K.carrier] - 1} := by
  ext n
  constructor
  · exact fun hn => ⟨hn.1, dvd_residueCard_sub_one_of_mem_spec K hn⟩
  · rintro ⟨hn0, hdvd⟩
    exact ⟨hn0, not_dvd_of_dvd_residueCard_sub_one K hdvd,
      contHomCard_absGal_of_dvd K (not_dvd_of_dvd_residueCard_sub_one K hdvd) hn0 hdvd⟩

/-! ## 輸送 -/

/-- `S` は位相群の同型で不変(`contHomCard_congr` の言い換え)。 -/
theorem residueCardSpec_congr {K K' : PAdicLocalField p}
    (α : K.absGal ≃ₜ* K'.absGal) : residueCardSpec K = residueCardSpec K' := by
  ext n
  simp only [residueCardSpec, Set.mem_setOf_eq, contHomCard_congr α n]

/-- **★★★★★★★★★★★★★★★★★★★★★★★★剰余体の位数 `q` は `Γ_K` から
群論的に決まる**。

[pGC] Proposition 1.2 の「`q` の回復」を、経路 C(`Interface` を経由しない道)で
確定させたもの。`Γ_K ≃ₜ* Γ_{K'}` なる位相群の同型があれば `#𝓀[K] = #𝓀[K']`。

証明:`S` が同型で不変(`residueCardSpec_congr`)で、その最大元が `q − 1`
(`isGreatest_residueCardSpec`)であり、最大元は一意(`IsGreatest.unique`)。
最後に `q, q' ≥ 2` から `q − 1 = q' − 1 ⟹ q = q'`。

★この定理は `inertia` にも `Interface` の `ResidueCardinality` にも依存しない。
依存は `contHomCard`(位相群だけで書ける不変量)と、その値を計算する
`Found/PGC` の下界・上界だけである。 -/
theorem residueCard_eq_of_absGal_equiv {K K' : PAdicLocalField p}
    (α : K.absGal ≃ₜ* K'.absGal) :
    Nat.card 𝓀[K.carrier] = Nat.card 𝓀[K'.carrier] := by
  have hgt : Nat.card 𝓀[K.carrier] - 1 = Nat.card 𝓀[K'.carrier] - 1 :=
    (isGreatest_residueCardSpec K).unique
      (by rw [residueCardSpec_congr α]; exact isGreatest_residueCardSpec K')
  have h1 : 1 < Fintype.card 𝓀[K.carrier] := Fintype.one_lt_card
  have h2 : 1 < Fintype.card 𝓀[K'.carrier] := Fintype.one_lt_card
  rw [Nat.card_eq_fintype_card, Nat.card_eq_fintype_card] at hgt ⊢
  omega

end ABC3.Found.PGC
