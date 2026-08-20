---
name: frdi-vloc-refuted
description: 𝒟 は 𝒞 の base-isomorphism による局所化ではない。反証は Example 4.3 で機械検証済み。原典は粗化 R_i を使う。
metadata:
  type: project
---

★2026-08-19、`Theorem 3.4, (v)` の `Ψ_Base` を組むために採った
**「`𝒟` は `𝒞` を base-isomorphism で局所化したもの」という枠は偽**である。
反証: `ABC3/Check/FrdI/VLocFalse.lean` の `ex43_not_isLocalization`(機械検証済み)。

## ★反例
`Example 4.3` の `𝒞`(在庫の `ex43P` / `ex43_core`)。
`𝒟 = Discrete PUnit` なので射はすべて base-iso。対象 `0 ∈ ℚ` の自己射は `ℕ≥1` 全体。
2 進付値 `ℕ≥1 →* Multiplicative ℤ` が `𝒞` から 1 対象亜群への関手を与え、
すべての射を反転するのに次数 2 の自己射を `1` へ送らない。

## ★★読み違えの型 —— **「充満」を「同値」と読んだ**
原典 (物理 p.68) の `φ_𝒟 = Base(ψ) ◦ Base(γ) ◦ Base(α)^{-1}` は
**充満性**しか言っていない。忠実性は別で、実際に**偽**である ——
`𝒪^×(A)` の元と Frobenius 次数は `Base` が潰すのに局所化には残る。

**Why:** 「分解の式が書ける」＝「局所化である」と即断した。
★局所化は**普遍性**の主張なので、分解(＝充満性)だけでは足りない。

**How to apply:**
1. 「X は Y の局所化」と読み替えたくなったら、**忠実性の側を必ず別に確かめる**。
   具体的には「**核に何が残るか**」——ここでは `𝒪^×` と Frobenius 次数。
2. ★★**在庫の具体例で反証を試す**のが速い。本プロジェクトは
   `Example43.lean` に `FrobenioidCore` を満たす具体圏を持っており、
   `𝒟` が自明なので**局所化系の主張はここで一撃で落ちる**。
   新しい枠を採ったら、まず `ex43P` に当てること。
3. 原典が付録(`Definition A.1` / `Proposition A.2`)を用意しているときは、
   **それを迂回する読み替えを疑う**。付録は飾りではない。

関連: [[frdi-split-nonisotropic-not-derivable]] [[no-wall-decompose-instead]]
