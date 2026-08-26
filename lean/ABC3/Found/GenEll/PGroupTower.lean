import ABC3.Meta.Claim
import Mathlib.GroupTheory.Sylow

/-!
# [GenEll] Proposition 1.7 の "elementary claim" の段 3 —— **`p`-群の塔**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.10。

原文 (GenEll p.10):
> Σ” of log-condE, log-condD is ≈0 [cf. Remark 1.5.1], while [again by the elementary

## ★★★★段 3 が言っていること

> Moreover, since `Gal(L/K)` is a [necessarily **solvable**!] `p`-group of order `≤ d`,
> it suffices to consider the case where `[L : K] = p`.

★これは Galois 対応で**体の塔** `K = K₀ ⊆ K₁ ⊆ … ⊆ K_s = L`(各段の次数が `p`)に
落とすということであり、群の側では**指数 `p` ずつ上がる部分群の鎖**である。

★★`different` の側は塔で継げる(`DifferentKummer.lean` の
`pow_mem_differentIdeal_tower`)ので、あとはこの鎖があればよい。

## ★★★在庫を引いた結果 —— mathlib の 1 本で足りる

`Sylow.exists_subgroup_card_pow_prime_le`:

> `p^m ∣ |G|` かつ `|H| = p^n` かつ `n ≤ m` なら、`H ≤ K` で `|K| = p^m` なる `K` がある

★これを `m ≝ n+1` で使えば**1 段上がる**。鎖はそれを繰り返すだけである。
★★**可解性を経由する必要はなかった**——原文が「necessarily solvable!」と書いたのは
読者への注意であって、鎖の存在には Sylow の定理で足りる。
-/

namespace ABC3.Found.GenEll

/-! ## ★★★★指数 `p` で 1 段上がる -/

/-- ★★★★**位数 `p^n` の部分群は、位数 `p^{n+1}` の部分群に 1 段だけ持ち上がる**。

原文 (GenEll p.10):
> Σ” of log-condE, log-condD is ≈0 [cf. Remark 1.5.1], while [again by the elementary

★原文の段 3(`Gal(L/K)` は位数 `≤ d` の `p`-群なので `[L:K] = p` の場合に帰着する)は、
Galois 対応でこの鎖に対応する。 -/
theorem exists_subgroup_succ_le {G : Type*} [Group G] [Finite G] {p : ℕ} [Fact p.Prime] {k : ℕ}
    (hcard : Nat.card G = p ^ k) (H : Subgroup G) {n : ℕ} (hn : Nat.card H = p ^ n)
    (hlt : n < k) :
    ∃ K : Subgroup G, Nat.card K = p ^ (n + 1) ∧ H ≤ K := by
  refine Sylow.exists_subgroup_card_pow_prime_le p ?_ H hn (by omega)
  rw [hcard]
  exact pow_dvd_pow p (by omega)

/-- ★★★**任意の中間の位数まで一気に上がる形**(上の繰り返し)。

★`m` 段上げたいときはこちらを直接使う。 -/
theorem exists_subgroup_card_pow_le {G : Type*} [Group G] [Finite G] {p : ℕ} [Fact p.Prime]
    {k : ℕ} (hcard : Nat.card G = p ^ k) (H : Subgroup G) {n : ℕ} (hn : Nat.card H = p ^ n)
    {m : ℕ} (hnm : n ≤ m) (hmk : m ≤ k) :
    ∃ K : Subgroup G, Nat.card K = p ^ m ∧ H ≤ K :=
  Sylow.exists_subgroup_card_pow_prime_le p (by rw [hcard]; exact pow_dvd_pow p hmk) H hn hnm

/-! ## ★出典の紐付け(`.src`) -/

def exists_subgroup_succ_le.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 10,
    item := "Proposition 1.7, (i) の elementary claim(段 3——p-群を指数 p の鎖に分ける)",
    sectionId := "genell-prop-1-7" }

def exists_subgroup_card_pow_le.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 10,
    item := "Proposition 1.7, (i) の elementary claim(段 3——鎖を一気に上げる形)",
    sectionId := "genell-prop-1-7" }

end ABC3.Found.GenEll
