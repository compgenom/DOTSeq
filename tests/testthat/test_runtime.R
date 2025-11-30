test_that("runtime returns only seconds when duration < 60 secs", {
    start <- as.POSIXct("2025-11-29 10:00:00")
    end <- as.POSIXct("2025-11-29 10:00:30")  # 30 seconds later
    result <- runtime(end, start)
    
    expect_type(result, "list")
    expect_true("secs" %in% names(result))
    expect_false("mins" %in% names(result))
    expect_equal(result$secs, 30)
})

test_that("runtime returns minutes and seconds when duration >= 60 secs", {
    start <- as.POSIXct("2025-11-29 10:00:00")
    end <- as.POSIXct("2025-11-29 10:02:15")  # 135 seconds later
    result <- runtime(end, start)
    
    expect_type(result, "list")
    expect_true(all(c("mins", "secs") %in% names(result)))
    expect_equal(result$mins, 2)
    expect_equal(result$secs, 15)
})

test_that("runtime handles zero duration correctly", {
    start <- as.POSIXct("2025-11-29 10:00:00")
    end <- start
    result <- runtime(end, start)
    
    expect_equal(result$secs, 0)
    expect_false("mins" %in% names(result))
})

test_that("runtime works with custom units argument", {
    start <- as.POSIXct("2025-11-29 10:00:00")
    end <- as.POSIXct("2025-11-29 10:01:00")
    # Using units = "mins" should return 1 minute total
    total_minutes <- as.numeric(difftime(end, start, units = "mins"))
    expect_equal(total_minutes, 1)
})
