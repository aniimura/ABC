/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.TensorPowCoherence
import ABC3.Meta.Claim

/-!
# ★★★★★★★★共通次数 `L` の有限個の切断で覆える —— **段 E2 が閉じた**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.6。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

## ★★★★★★★★これで段 E2 が閉じた

    §9-817  有限被覆（準コンパクト性）
    §9-818  `X_{s^{⊗k}} = X_s`
    §9-819  `M^{⊗(n·k)} ≅ (M^{⊗n})^{⊗k}`（コヒーレンス）
    ★本ファイル  `lcm` の帳簿 —— **共通次数 `L` の切断で覆える**

## ★★★★★機構 —— `lcm` ではなく**積**でよい

`L = ∏_{c ∈ T} c.n` と取れば、各 `c.n ∣ L` は `Finset.dvd_prod_of_mem` で出る。
★`lcm` を取る必要はない——**倍数でありさえすればよい**からである。
★★`L > 0` は `Finset.prod_pos`（各 `c.n > 0`）。

★★★各 `c` について `k_c = L / c.n` とし、`secPowAligned`（本ファイル）が
`s_c^{⊗k_c}` を `M^{⊗L}` の中へ運ぶ。非消失軌跡は変わらない（`§9-818` ＋ `§9-819` ＋
`nonVanishing_of_iso`）。

## ★測定の記録

`c.n ∣ L` から `L = c.n * k` を取り出して **`subst`** すれば、
`M^{⊗L}` と `M^{⊗(c.n·k)}` の型の食い違いは消える。
★`▸` で運ぶ必要はなかった（2026-08-28 実測）。

## ★残っている段（明示）

★★段 E2 は閉じた。残るのは **E1d（貼り合わせの重なり一致）**と
**E3（immersion 性）**、および E0 の等号 1 本・C2c 前半・F2c である。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite ABC3.Found.Arakelov

variable {X : Scheme.{0}}

/-! ## ★★★★★(1) 次数を揃えた切断 -/

/-- ★★★★★**次数を `n·k` に揃えた切断** `s^{⊗k}`。

★`secPow` の行き先は `(M^{⊗n})^{⊗k}` なので、`§9-819` のコヒーレンスで
`M^{⊗(n·k)}` へ運ぶ。 -/
noncomputable def secPowAligned (M : X.PresheafOfModules) (n k : ℕ)
    (s : ((presheafTensorPow M n).obj (op ⊤) : Type)) :
    ((presheafTensorPow M (n * k)).obj (op ⊤) : Type) :=
  ((presheafTensorPowMul M n k).symm).hom.app (op ⊤) (secPow (presheafTensorPow M n) s k)

/-- ★★★★★★**次数を揃えても非消失軌跡は変わらない**。 -/
theorem nonVanishing_secPowAligned (M : X.PresheafOfModules) (hM : IsLocallyTrivial X M)
    (n k : ℕ) (s : ((presheafTensorPow M n).obj (op ⊤) : Type)) :
    nonVanishing (presheafTensorPow M (n * (k + 1))) (secPowAligned M n (k + 1) s)
      = nonVanishing (presheafTensorPow M n) s := by
  rw [secPowAligned,
    nonVanishing_of_iso _ _ ((presheafTensorPowMul M n (k + 1)).symm),
    nonVanishing_secPow (presheafTensorPow M n)
      (isLocallyTrivial_presheafTensorPow hM n) s k]

/-- ★★★★★★**任意の共通次数 `L`（各 `c.n` の倍数）へ切断を持ち上げられる**。

★`L = c.n * k` を `subst` すれば型の食い違いは消える。 -/
theorem exists_section_of_degree (M : X.PresheafOfModules) (hM : IsLocallyTrivial X M)
    (c : AmpleChart M) (L : ℕ) (hdvd : c.n ∣ L) (hL : 0 < L) :
    ∃ t : ((presheafTensorPow M L).obj (op ⊤) : Type),
      nonVanishing (presheafTensorPow M L) t = c.open' := by
  obtain ⟨k, hk⟩ := hdvd
  subst hk
  have hkpos : 0 < k := by
    rcases Nat.eq_zero_or_pos k with h | h
    · rw [h, Nat.mul_zero] at hL; exact absurd hL (lt_irrefl 0)
    · exact h
  obtain ⟨k', rfl⟩ : ∃ k', k = k' + 1 := ⟨k - 1, by omega⟩
  exact ⟨secPowAligned M c.n (k' + 1) c.s, nonVanishing_secPowAligned M hM c.n k' c.s⟩

/-! ## ★★★★★★★★(2) 段 E2 の結論 -/

open scoped Classical in
set_option maxHeartbeats 1000000 in
/-- ★★★★★★★★**ample なら、ある共通次数 `L` の有限個の切断で `X` を覆える**。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

★★これが射 `X ⟶ ℙᴺ_R` を作るための**入力**である
——同じ束 `M^{⊗L}` の有限個の切断が `X` を覆う。
★★★`L` は `lcm` ではなく**積** `∏ c.n` で十分である（倍数でありさえすればよい）。 -/
theorem exists_common_degree_cover (M : X.PresheafOfModules) (hM : IsLocallyTrivial X M)
    (h : IsAmple M) :
    ∃ L : ℕ, 0 < L ∧ ∃ S : Set ((presheafTensorPow M L).obj (op ⊤) : Type), S.Finite ∧
      (⨆ u ∈ S, (nonVanishing (presheafTensorPow M L) u : Set X)) = Set.univ ∧
      ∀ u ∈ S, IsAffineOpen (nonVanishing (presheafTensorPow M L) u) := by
  obtain ⟨T, hTfin, hTcov⟩ := exists_finite_cover_of_isAmple M h
  set F := hTfin.toFinset with hF
  set L := ∏ c ∈ F, c.n with hLdef
  have hLpos : 0 < L := Finset.prod_pos (fun c _ => c.hn)
  choose t ht using fun (c : F) => exists_section_of_degree M hM c.1 L
    (Finset.dvd_prod_of_mem _ c.2) hLpos
  refine ⟨L, hLpos, Set.range t, Set.finite_range t, ?_, ?_⟩
  · refine Set.eq_univ_of_univ_subset ?_
    rw [← hTcov]
    intro x hx
    simp only [Set.iSup_eq_iUnion, Set.mem_iUnion, exists_prop] at hx ⊢
    obtain ⟨c, hcT, hxc⟩ := hx
    have hcF : c ∈ F := hTfin.mem_toFinset.2 hcT
    refine ⟨t ⟨c, hcF⟩, ⟨⟨c, hcF⟩, rfl⟩, ?_⟩
    rw [ht ⟨c, hcF⟩]
    exact hxc
  · rintro u ⟨c, rfl⟩
    rw [ht c]
    exact c.1.isAffineOpen_open'

/-! ## ★出典の紐付け(`.src`) -/

def secPowAligned.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(次数を n·k に揃えた切断)",
    sectionId := "genell-prop-1-4" }

def exists_section_of_degree.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(共通次数へ切断を持ち上げられる)",
    sectionId := "genell-prop-1-4" }

def exists_common_degree_cover.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(ample なら共通次数 L の有限個の切断で X を覆える)",
    sectionId := "genell-prop-1-4" }

def exists_common_degree_cover.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "exists_finite_cover_of_isAmple(有限被覆、§9-817)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_finite_cover_of_isAmple") 6,
    .citation "[ABC3]" "nonVanishing_secPow(X_{s^{⊗k}} = X_s、§9-818)"
      (.inProject "ABC3" "ABC3.Found.GenEll.nonVanishing_secPow") 6,
    .citation "[ABC3]" "presheafTensorPowMul(コヒーレンス、§9-819)"
      (.inProject "ABC3" "ABC3.Found.GenEll.presheafTensorPowMul") 6,
    .citation "[mathlib]" "Finset.dvd_prod_of_mem / Finset.prod_pos"
      (.inMathlib "Finset.dvd_prod_of_mem") 6,
    .implicitStep
      ("★L は lcm ではなく**積** ∏ c.n で十分である——倍数でありさえすればよい。" ++
       "★★c.n ∣ L から L = c.n * k を取り出して subst すれば型の食い違いは消える" ++
       "(▸ で運ぶ必要はなかった、2026-08-28 実測)") 6,
    .implicitStep
      ("★★★段 E2 は閉じた。残るのは E1d(貼り合わせの重なり一致)と E3(immersion 性)、" ++
       "および E0 の等号 1 本・C2c 前半・F2c である") 6 ]

end ABC3.Found.GenEll
